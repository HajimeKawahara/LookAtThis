import pandas as pd
import numpy as np
import typetomass as ttm
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_sun
import matplotlib.pyplot as plt
import argparse

def is_float_expression(s):
    try:
        val = float(s)
        return True
    except ValueError:
        return False
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Look at this from SPOTS!')
    parser.add_argument('-t', nargs=2, default=[100.0,2000.0],help='separation min and max (mas)',type=float)
    parser.add_argument('-a', nargs=1, default=[40.0],help='maximum elevation [deg]',type=float)
    parser.add_argument('-p', nargs=3, default=[19.828611,155.48055,4139.],help='latitude/longitude/height(m) of the observatory. Default is Maunakea',type=float)
    parser.add_argument('-d', nargs=1, default=["2019-7-15"],help='observation date',type=float)

    args = parser.parse_args()    


    
    thetamin=args.t[0]#100.0 #mas
    thetamax=args.t[1]#2000.0 #mas
    maxalt = args.a[0]#30.0 #maximum altitude
    lat_subaru=args.p[0]#19.0 + 49/60.+ 43/3600.0
    lon_subaru=args.p[1]#155.0 + 28/60.0 + 50/3600.0 
    height_subaru = args.p[2]#4139.0
    utcoffset = - 10.0*u.hour
    midlocal=Time(args.d)+(1*u.d)
    #midlocal=Time('2019-7-16 00:00:00')
    midnight = midlocal + utcoffset
    location_subaru=EarthLocation(lat=lat_subaru*u.deg, lon=lon_subaru*u.deg, height=height_subaru*u.m)

    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_obs = midnight + delta_midnight
    frame=AltAz(obstime=midnight+delta_midnight,location=location_subaru)
    sunaltazs = get_sun(midnight+delta_midnight).transform_to(frame)
    nightmask=(sunaltazs.alt < -0*u.deg)

    
    dat=pd.read_csv("../database/spots/spots.tsv",delimiter="|")

    fig=plt.figure(figsize=(15,7))
    ax=fig.add_subplot(121)
    ic=0
    lsarr=["solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted"]
    for i,name in enumerate(dat["Star"]):
        d=float(dat["Dist"][i])
        sp=dat["SpType"][i]

        mass1=float(dat["MA"][i])
        try:
            mass2=float(dat["MB"][i])
            totalmass=mass1+mass2
        except:
            totalmass=mass1
            
        acrit=float(dat["acrit"][i])
        if is_float_expression(dat["rho"][i]):
            theta=float(dat["rho"][i])*1000
        else:
            try:
                P=float(dat["Per"][i])
                if dat["x_Per"][i] == "yr":
                    P=P*365
                    
                a=((P/365.0)**2*totalmass)**(1.0/3.0)
                theta=a/d*1000
            except:
                theta=-1.0
                

        
        #print(name,c.to_string('hmsdms'))
        if theta > thetamin and theta < thetamax:
            print(name,sp,mass1,mass2,"theta=",theta)
            ra=dat["RAJ2000"][i]
            dec=dat["DEJ2000"][i]
            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

            altitude = c.transform_to(frame)
            if np.max(altitude.alt[nightmask]) > maxalt*u.deg:
                iic=np.mod(ic,7)
                lab=name+", "+sp+", mA="+str(mass1)+", mB="+str(mass2)+", $\\theta$= "+str(round(theta,1))+" mas,"+" H="+str(round(float(dat["Hmag"][i]),1))
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
    plt.title("SPOTS")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
    plt.show()
