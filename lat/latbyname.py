import pandas as pd
import numpy as np
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
    parser = argparse.ArgumentParser(description='Just show elevation for names given by a user')
    parser.add_argument('-n', nargs="+", help='name',type=str)
    parser.add_argument('-a', nargs=1, default=[40.0],help='maximum elevation [deg]',type=float)
    parser.add_argument('-p', nargs=3, default=[19.828611,204.51945,4139.],help='latitude/longitude/height(m) of the observatory. Default is Maunakea',type=float)
    parser.add_argument('-d', nargs=1, default=["2019-7-15"],help='observation date',type=str)

    args = parser.parse_args()        
    maxalt = args.a[0]#30.0 #maximum altitude
    lat_subaru=args.p[0]#19.0 + 49/60.+ 43/3600.0
    lon_subaru=args.p[1]#155.0 + 28/60.0 + 50/3600.0 
    height_subaru = args.p[2]#4139.0
    utcoffset = 10.0*u.hour
    midlocal=Time(args.d)+(1*u.d)
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
#    nightmask=(sunaltazs.alt < -0*u.deg)
#    print(len(delta_midnight[nightmask]))

    nightmask=(sunaltazs.alt < -0*u.deg)*maskobs

    fig=plt.figure(figsize=(15,7))
    ax=fig.add_subplot(121)
    ic=0
    lsarr=["solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted"]
    for i,name in enumerate(args.n):
        if True:
            Simbad.SIMBAD_URL = "http://simbad.u-strasbg.fr/simbad/sim-script"

            Simbad.add_votable_fields("sp","flux(V)","flux(R)","flux(J)","flux(H)","flux(K)")                
            result_table = Simbad.query_object(name)
            print(result_table)
            namex=result_table["MAIN_ID"][0].decode('utf-8')
            ra=result_table["RA"][0]
            dec=result_table["DEC"][0]
            V=result_table["FLUX_V"].item()
            R=result_table["FLUX_R"].item()        
            J=result_table["FLUX_J"].item()
            H=result_table["FLUX_H"].item()
            K=result_table["FLUX_K"].item()            
            
            sptype=result_table["SP_TYPE"][0].decode('utf-8')
            print("==>",name,ra,dec,sptype)


            c = SkyCoord(ra+" "+dec, unit=(u.hourangle, u.deg))
            #            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

            altitude = c.transform_to(frame)
            iic=np.mod(ic,7)
            lab=name+" ("+namex+"), "+(sptype)+" V="+str(round(V,1))+" R="+str(round(R,1))+" J="+str(round(J,1))+" H="+str(round(H,1))+" K="+str(round(K,1))
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
            plt.xlabel('Hours from Midnight = '+midlocal.iso[0])
            plt.ylabel('Altitude [deg]')
            plt.title("")
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
    plt.show()
