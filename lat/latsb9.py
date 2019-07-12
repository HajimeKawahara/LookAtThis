import pandas as pd
import numpy as np
import typetomass as ttm
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_sun
import matplotlib.pyplot as plt
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Look at this from SPOTS!')
    parser.add_argument('-t', nargs=2, default=[100.0,2000.0],help='separation min and max (mas)',type=float)
    parser.add_argument('-a', nargs=1, default=[40.0],help='maximum elevation [deg]',type=float)
    parser.add_argument('-p', nargs=3, default=[19.828611,155.48055,4139.],help='latitude/longitude/height(m) of the observatory. Default is Maunakea',type=float)
    parser.add_argument('-d', nargs=1, default=["2019-7-15"],help='observation date',type=str)
    parser.add_argument('-m', nargs=1, default=["second night"],help='observing mode: full night (default), first night, second night',type=str)

    args = parser.parse_args()        
    thetamin=args.t[0]#100.0 #mas
    thetamax=args.t[1]#2000.0 #mas
    maxalt = args.a[0]#30.0 #maximum altitude
    lat_subaru=args.p[0]#19.0 + 49/60.+ 43/3600.0
    lon_subaru=args.p[1]#155.0 + 28/60.0 + 50/3600.0 
    height_subaru = args.p[2]#4139.0
    utcoffset = - 10.0*u.hour
    midlocal=Time(args.d)+(1*u.d)
    midnight = midlocal + utcoffset

    obsmode=args.m[0]
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


    
    dat=pd.read_csv("../database/sb9/sb9.csv",comment="#")
    dbV=pd.read_csv("../database/sp/spV.txt",delimiter=",")
#    print(dbV["mass"][dbV["Sp"]=="O8"].values[0])

    fig=plt.figure(figsize=(15,7))
    ax=fig.add_subplot(121)
    ic=0
    lsarr=["solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted"]
    for i,sp1 in enumerate(dat["Sp1"]):
        sysnum=dat["System number"][i]
        mass1=ttm.type2mass(sp1,dbV,unknownmass=0.0)
        sp2=dat["Sp2"][i]
        mass2=ttm.type2mass(sp2,dbV,unknownmass=0.0)
        totalmass=mass1+mass2
        P=dat["Period (d)"][i]
        a=((P/365.0)**2*totalmass)**(1.0/3.0)
        plx_simbad=dat["Simbad plx"][i]
        plx_gaia=dat["GAIA plx"][i]
        theta_gaia=a*(plx_gaia)
        theta_simbad=a*(plx_simbad)
        K1=dat["K1 (km/s)"][i]
        K2=dat["K2 (km/s)"][i]

        if sp2 != sp2:
            sp2="N/A"
        #print(sp1,mass1,sp2,mass2)
        #print("Mtot=",totalmass,"P=",P,"(d) a=",a,"au")
        #print("d=",1/(plx_simbad/1000),"pc")
        #print("theta gaia=",theta_gaia,"mas","theta simbad=",theta_simbad,"mas")
        
        if theta_simbad > thetamin and theta_simbad < thetamax and totalmass > 0.0 and mass2 == 0.0:
            name=dat["name"][i]
            radec=dat["radec"][i]
            c=SkyCoord.from_name("J"+radec, parse=True)
            altitude = c.transform_to(frame)
            if np.max(altitude.alt[nightmask]) > maxalt*u.deg:
                print(name,sp1,mass1,sp2,mass2,"theta simbad=",theta_simbad,"mas",K1,K2)
                iic=np.mod(ic,8)
                print(sp1,sp2,mass1,mass2)
                lab=name+", "+sp1+","+sp2+", mA="+str(mass1)+", mB="+str(mass2)+", $\\theta$= "+str(round(theta_simbad,1))+" mas,"+" H="+str(round(float(dat["H"][i]),1))+",KA="+str(round(K1,0))+" km/s, KB="+str(round(K2,0))+" km/s"
                plt.plot(delta_midnight,altitude.alt,label=lab,color="C"+str(iic),ls=lsarr[int(ic/8)])
                ic=ic+1

            #            print("Mtot=",totalmass,"P=",P,"(d) a=",a,"au")
#            print("d=",1/(plx_simbad/1000),"pc")
#            print("theta gaia=",theta_gaia,"mas","theta simbad=",theta_simbad,"mas")
            

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
    plt.title("SB9")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
    plt.show()

