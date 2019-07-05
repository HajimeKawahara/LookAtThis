import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

dat=pd.read_csv("../sb9/Main.dta",delimiter="|",dtype={"System Number":"int","1900.0 coordinates":"str","2000.0 coordinates":"str","Component":"str","Magnitude of component 1":"float","Filter component 1":"str","Magnitude of component 2":"float","Filter component 2":"str","Spectral type component 1":"str","Spectral type component 2":"str"})

dat2=pd.read_csv("../sb9/Orbits.dta",delimiter="|",dtype={"System number":"int","Orbit number for that system":"int","Period (d)":"float","error on P (d)":"float","Periastron time (JD-2400000)":"float","error on Periastron time":"float","Flag on periastron time":"str","eccentricity":"float","error on eccentricity":"float","argument of periastron (deg)":"float","error on omega":"float","K1 (km/s)":"float","error on K1 (km/s)":"float","K2 (km/s)":"float","error on K2 (km/s)":"float","systemic velocity (km/s)":"float","error on V0 (km/s)":"float","rms RV1 (km}":"float","rms RV2 (km/s)":"float","RV1":"float","RV2":"float","Grade":"float","Bibcode":"str","Contributor":"str","Accessibility":"str","Reference adopted for the times (JD or MJD)":"str"},na_values = "NaN")



dat4=pd.read_csv("../sb9/sb9_position_clean.txt",delimiter="|",dtype={"System Number":"int","ra (degree)":"float","dec (degree)":"float","Simbad plx":"float","GAIA plx":"float","V":"float","R":"float","J":"float","H":"float","K":"float"})
dat4=dat4.set_index("System Number")

newd=dat2.set_index("System number")
newd["Sp1"]="" 
newd["Sp2"]="" 
newd["radec"]= ""
newd["Simbad plx"]= np.nan
newd["GAIA plx"]= np.nan
newd["V"]= np.nan
newd["R"]= np.nan
newd["J"]= np.nan
newd["H"]= np.nan
newd["K"]= np.nan
newd["name"]=""

dat3=pd.read_csv("../sb9/Alias.dta",delimiter="|",dtype={"System number":"int","Catalog name":"str","ID in that catalog":"str"})

for i,sysi in enumerate(dat["System Number"]):
    try:
        #|Catalog name|ID in that catalog
        mask=dat3["System number"]==sysi
        ID=dat3["ID in that catalog"][mask].values
        for jj,tag in enumerate(dat3["Catalog name"][mask].values):
            if jj == 0:
                name=tag+" "+ID[jj]
            if tag == "HD":
                name=tag+" "+ID[jj]
        newd.loc[sysi,"name"]=name
        #print(i,sysi)
        newd.loc[sysi,"Sp1"]=dat["Spectral type component 1"][i]
        newd.loc[sysi,"Sp2"]=dat["Spectral type component 2"][i]
        
        newd.loc[sysi,"radec"]=dat["2000.0 coordinates"][i]
        newd.loc[sysi,"Simbad plx"]=dat4.loc[sysi,"Simbad plx"]
        newd.loc[sysi,"GAIA plx"]=dat4.loc[sysi,"GAIA plx"]
        newd.loc[sysi,"V"]=dat4.loc[sysi,"V"]
        newd.loc[sysi,"R"]=dat4.loc[sysi,"R"]
        newd.loc[sysi,"J"]=dat4.loc[sysi,"J"]
        newd.loc[sysi,"H"]=dat4.loc[sysi,"H"]
        newd.loc[sysi,"K"]=dat4.loc[sysi,"K"]
    except:
        print("No data for ",sysi)
    

newd.to_csv("sb9.csv")
