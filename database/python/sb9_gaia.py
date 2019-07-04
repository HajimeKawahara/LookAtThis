import pandas as pd
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

dat3=pd.read_csv("../database/sb9/Alias.dta",delimiter="|",dtype={"System number":"int","Catalog name":"str","ID in that catalog":"str"})


for i in range(0,len(dat3["Catalog name"])):
    if dat3["Catalog name"][i]!="GAIADR2":
        print(dat3["Catalog name"][i]+" "+dat3["ID in that catalog"][i])
