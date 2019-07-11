import pandas as pd
import numpy as np
import typetomass as ttm
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_sun
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
import sys
#Seq|Star|RAJ2000|DEJ2000|Dist|n_Dist|Age|Hmag|SpType|MA|l_MB|MB|l_rho|rho|l_Per|Per|x_Per|l_Ecc|Ecc|l_acrit|acrit|Ref|SimbadName

def is_float_expression(s):
    try:
        val = float(s)
        return True
    except ValueError:
        return False
    
if __name__ == "__main__":

    obsmode="second night"

#    low contrast
#    contmax=0.1
#    contmin=0.05

    
    #mid contrast
    contmax=0.05
    contmin=0.03

    #high contrast
#    contmax=0.03
#    contmin=0.001

    
    thetamin=100.0 #mas
    thetamax=2000.0 #mas
    maxalt = 40.0 #maximum altitude
    lat_subaru=19.0 + 49/60.+ 43/3600.0
    lon_subaru=155.0 + 28/60.0 + 50/3600.0 
    height_subaru = 4139.0
    utcoffset = - 10.0*u.hour
    midlocal=Time('2019-7-16 00:00:00')
    midnight = midlocal + utcoffset

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
#    nightmask=(sunaltazs.alt < -0*u.deg)
#    print(len(delta_midnight[nightmask]))

    nightmask=(sunaltazs.alt < -0*u.deg)*maskobs
    
    dat=pd.read_csv("../database/siwiyn/siwiyn.csv",delimiter=",",comment="#")
 #,_RAJ2000,_DEJ2000,WDS,Name,HD,HIP,Date,PA,Sep,Dmag,Wave,FWHM,f_FWHM,SimbadName,_RA,_DE,Simbad plx,GAIA plx,V,R,J,H,K

    fig=plt.figure(figsize=(15,7))
    ax=fig.add_subplot(121)
    ic=0
    lsarr=["solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted"]
    for i,name in enumerate(dat["Name"]):
        theta=float(dat["Sep"][i])*1000.0 #[mas]
        try: 
            contrast = 1.0/10**(0.4*float(dat["Dmag"][i]))
        except:
            contrast=0.0
            #print(name,c.to_string('hmsdms'))
            #        if True:
        if theta > thetamin and theta < thetamax and contrast < contmax and contrast > contmin:
            ra=dat["_RAJ2000"][i]
            dec=dat["_DEJ2000"][i]
            c = SkyCoord(ra+" "+dec, unit=(u.hourangle, u.deg))
            #            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))



            altitude = c.transform_to(frame)
            if np.max(altitude.alt[nightmask]) > maxalt*u.deg:
                print(name,contrast)

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
                    except:
                        print("==")
                
                try:
                    Simbad.SIMBAD_URL = "http://simbad.u-strasbg.fr/simbad/sim-script"
                    Simbad.add_votable_fields("sp")                
                    result_table = Simbad.query_object(name)
                    sptype=result_table["SP_TYPE"][0].decode('utf-8')
                except:
                    sptype="N/A"
                iic=np.mod(ic,7)
                #                lab=name+", "+sp+", mA="+str(mass1)+", mB="+str(mass2)+", $\\theta$= "+str(round(theta,1))+" mas,"+" H="+str(round(float(dat["Hmag"][i]),1))
                lab=name+", "+(sptype)+", contrast="+str(round(contrast,3))+", $\\theta$= "+str(round(theta,1))+" mas,"+" V="+str(round(V,1))+" R="+str(round(R,1))+" H="+str(round(H,1))
                plt.plot(delta_midnight,altitude.alt,label=lab,color="C"+str(iic),ls=lsarr[int(ic/7)])
                ic=ic+1
                
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                 sunaltazs.alt < -0*u.deg, color='0.5', zorder=0)
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                 sunaltazs.alt < -18*u.deg, color='k', zorder=0)
    plt.ylim(0.,90)
    plt.xlim(-12, 12)
    plt.xticks(np.arange(13)*2 -12)
    plt.ylim(10, 90)
    plt.axhline(30.0,color="gray",lw=1)
    plt.xlabel('Hours from Midnight = '+midlocal.iso)
    plt.ylabel('Altitude [deg]')
    plt.title("Speckle imaging (WIYN)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
    plt.show()
