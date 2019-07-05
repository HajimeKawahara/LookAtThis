import pandas as pd
import numpy as np
import typetomass as ttm
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_sun
import matplotlib.pyplot as plt

if __name__ == "__main__":

    thetamin=100.0 #mas
    thetamax=1000.0 #mas
    maxalt = 30.0 #maximum altitude
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
    nightmask=(sunaltazs.alt < -0*u.deg)


    
    dat=pd.read_csv("../database/sb9/sb9.csv")
    dbV=pd.read_csv("../database/sp/spV.txt",delimiter=",")
#    print(dbV["mass"][dbV["Sp"]=="O8"].values[0])
    for i,sp1 in enumerate(dat["Sp1"]):
        mass1=ttm.type2mass(sp1,dbV,unknownmass=1.0)
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

        #print(sp1,mass1,sp2,mass2)
        #print("Mtot=",totalmass,"P=",P,"(d) a=",a,"au")
        #print("d=",1/(plx_simbad/1000),"pc")
        #print("theta gaia=",theta_gaia,"mas","theta simbad=",theta_simbad,"mas")
        
        if theta_simbad > thetamin and theta < thetamax:
            print(sp1,mass1,sp2,mass2,"theta simbad=",theta_simbad,"mas",K1,K2)
#            print("Mtot=",totalmass,"P=",P,"(d) a=",a,"au")
#            print("d=",1/(plx_simbad/1000),"pc")
#            print("theta gaia=",theta_gaia,"mas","theta simbad=",theta_simbad,"mas")
            
