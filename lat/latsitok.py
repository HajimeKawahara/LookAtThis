import pandas as pd
import numpy as np
import typetomass as ttm
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_sun
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
import argparse
import sys
#Seq|Star|RAJ2000|DEJ2000|Dist|n_Dist|Age|Hmag|SpType|MA|l_MB|MB|l_rho|rho|l_Per|Per|x_Per|l_Ecc|Ecc|l_acrit|acrit|Ref|SimbadName

def is_float_expression(s):
    try:
        val = float(s)
        return True
    except ValueError:
        return False
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Look at this from Speckel Images of Tokovinin19/WIYN!')
    parser.add_argument('-t', nargs=2, default=[100.0,2000.0],help='separation min and max (mas)',type=float)
    parser.add_argument('-a', nargs=1, default=[40.0],help='maximum elevation [deg]',type=float)
    parser.add_argument('-p', nargs=3, default=[19.828611,155.48055,4139.],help='latitude/longitude/height(m) of the observatory. Default is Maunakea',type=float)
    parser.add_argument('-d', nargs=1, default=["2019-9-07"],help='observation date',type=str)
    parser.add_argument('-m', nargs=1, default=["full night"],help='observing mode: full night (default), first night, second night',type=str)
#    parser.add_argument('-c', nargs=2, default=[0.03,0.1],help='contrast min max',type=float)
    parser.add_argument('-i', nargs=1, default=["info.txt"],help='info file',type=str)

    args = parser.parse_args()
    infofile=args.i[0]
    thetamin=args.t[0]#100.0 #mas
    thetamax=args.t[1]#2000.0 #mas
    maxalt = args.a[0]#30.0 #maximum altitude
    lat_subaru=args.p[0]#19.0 + 49/60.+ 43/3600.0
    lon_subaru=args.p[1]#155.0 + 28/60.0 + 50/3600.0

    today=args.d[0]
    timediff=Time(today)-Time("2000-1-01")
    yeardiff=(timediff.jd/365.242189)
    
    height_subaru = args.p[2]#4139.0
    utcoffset = - 11.0*u.hour
    midlocal=Time(today)+(1*u.d)
    midnight = midlocal + utcoffset

    obsmode=args.m[0]
    #mid contrast
#    contmax=args.c[1]
#    contmin=args.c[0]

    print(midlocal.iso)
    location_subaru=EarthLocation(lat=lat_subaru*u.deg, lon=lon_subaru*u.deg, height=height_subaru*u.m)
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_obs = midnight + delta_midnight    
    frame=AltAz(obstime=midnight+delta_midnight,location=location_subaru)
    sunaltazs = get_sun(midnight+delta_midnight).transform_to(frame)

    #### OBS mode ####
    maskobs=np.ones(len(delta_midnight),dtype=bool)
    ichangepoint=np.argmin(sunaltazs.alt.degree)
    if obsmode=="first night":
        print("first night")
        maskobs[ichangepoint:]=False
    elif obsmode=="second night":
        print("second night")
        maskobs[0:ichangepoint]=False
    else:
        print("full night")
#    nightmask=(sunaltazs.alt < -0*u.deg)
#    print(len(delta_midnight[nightmask]))

    nightmask=(sunaltazs.alt < -0*u.deg)*maskobs
    
    dat=pd.read_csv("../database/sitok/sitok_position.txt",delimiter="|",comment="#")

    fig=plt.figure(figsize=(15,7))
    ax=fig.add_subplot(121)
    ic=0
    lsarr=["solid","dashed","dotted","dashdot","solid","dashed","dotted","dashdot","solid","dashed","dotted","dashdot","solid","dashed","dotted","dashdot","solid","dashed","dotted","dashdot","solid","dashed","dotted","dashdot","solid","dashed","dotted","dashdot","solid","dashed","dotted"]
    ff=open(infofile,"a")
    for i,name in enumerate(dat["Name"]):
        theta=float(dat["sep"][i])*1000.0 #[mas]
            #print(name,c.to_string('hmsdms'))
            #        if True:
        if theta > thetamin and theta < thetamax:
            sp=dat["sp"][i]
            pradec=dat["radec"][i]
            c = SkyCoord("J"+pradec, unit=(u.hourangle, u.deg))
            #            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))


            altitude = c.transform_to(frame)
            if np.max(altitude.alt[nightmask]) > maxalt*u.deg:

                J=dat["J"][i]
                H=dat["H"][i]
                K=dat["K"][i]
                R=dat["R"][i]
                V=dat["V"][i]

                if float(H) != float(H):
                    try:
                        Simbad.SIMBAD_URL = "http://simbad.u-strasbg.fr/simbad/sim-script"
                        Simbad.add_votable_fields("flux(V)","flux(R)","flux(J)","flux(H)","flux(K)")

                        result_table = Simbad.query_object(name)
                        V=result_table["FLUX_V"].item()
                        R=result_table["FLUX_R"].item()        
                        J=result_table["FLUX_J"].item()
                        H=result_table["FLUX_H"].item()
                        K=result_table["FLUX_K"].item()
                        print(result_table)

                        sys.exit()
                    except:
                        print("==")
                
                sptype=sp
                print(R)
                raJ2000=dat["RA"][i]
                decJ2000=dat["DEC"][i]
                pmra=dat["PMRA"][i]
                pmdec=dat["PMDEC"][i]

                pradecJ2000=raJ2000+" "+decJ2000
                
#                cc=SkyCoord(pradecJ2000,unit=(u.hourangle, u.deg),pm_ra_cosdec=pmra * u.mas/u.yr,pm_dec=pmdec * u.mas/u.yr, obstime=Time(today))
#                print(cc.ra.degree+)
#http://www.sc.eso.org/~bdias/pycoffee/codes/20160908/pycoffee_coordinates.ipynb
                cc=SkyCoord(pradecJ2000,unit=(u.hourangle, u.deg))
                j2000=(cc.to_string('hmsdms'))
                ranow=(cc.ra+pmra*yeardiff/1000.0*u.arcsec/np.cos(cc.dec))
                decnow=(cc.dec+pmdec*yeardiff/1000.0*u.arcsec)
                print("YEAR DIFFERECE = ",yeardiff)
                print(pmra,"mas/yr")
                print(pmdec,"mas/yr")
                print(pmra*yeardiff/1000,"arcsec (RA)")
                print(pmdec*yeardiff/1000,"arcsec (DEC)")
                pradecnow=str(ranow)+" "+str(decnow)
                ccs=SkyCoord(ranow,decnow)
#                ccs=SkyCoord(pradecnow,unit=(u.deg, u.deg))
                nowpos=(ccs.to_string('hmsdms'))
                radecin=nowpos+","+j2000+","+str(pmra)+","+str(pmdec)
                if R != R:
                    ff.write(name+","+radecin+",V="+str(round(float(V),1))+"\n")
                    print(name,",",radecin,",V=",round(float(V),1))
                else:
                    ff.write(name+","+radecin+",R="+str(round(float(R),1))+"\n")
                    print(name,",",radecin,",R=",round(float(R),1))


                    
                iic=np.mod(ic,7)
                #                lab=name+", "+sp+", mA="+str(mass1)+", mB="+str(mass2)+", $\\theta$= "+str(round(theta,1))+" mas,"+" H="+str(round(float(dat["Hmag"][i]),1))
                lab=name+", "+(sptype)+", $\\theta$= "+str(round(theta,1))+" mas,"+" V="+str(round(float(V),1))+" R="+str(round(float(R),1))+" H="+str(round(float(H),1))
                plt.plot(delta_midnight,altitude.alt,label=lab,color="C"+str(iic),ls=lsarr[int(ic/7)])
                ic=ic+1
    ff.close()
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                 sunaltazs.alt < -0*u.deg, color='0.5', zorder=0)
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                 sunaltazs.alt < -18*u.deg, color='k', zorder=0)
    plt.ylim(0.,90)
    plt.xlim(-12, 12)
    plt.xticks(np.arange(13)*2 -12)
    plt.ylim(10, 90)
    plt.axhline(30.0,color="gray",lw=1)
    plt.xlabel('Hours from Midnight = '+midlocal.iso[0])
    plt.ylabel('Altitude [deg]')
    plt.title("Speckle imaging (Tokovinin19)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
    plt.show()
