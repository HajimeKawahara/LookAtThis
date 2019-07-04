import pandas as pd
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='get information from Simbad and GAIA TAP (DR1)')
parser.add_argument('-i', nargs=1, default=[0], help='start index', type=int)
parser.add_argument('-f', nargs=1, help='file', type=str)

args = parser.parse_args()    
istart = args.i[0]


dat=pd.read_csv("sb9/Main.dta",delimiter="|",dtype={"System Number":"int","1900.0 coordinates":"str","2000.0 coordinates":"str","Component":"str","Magnitude of component 1":"float","Filter component 1":"str","Magnitude of component 2":"float","Filter component 2":"str","Spectral type component 1":"str","Spectral type component 2":"str"})


#GAIA distance
#f=open("sb9_position.txt","a")
#f.write("System Number|ra (degree)|dec (degree)|Simbad plx|GAIA plx|V|R|J|H|K"+"\n")
#f.close()
if args.f:
    namelist=np.loadtxt(args.f[0],dtype=int)
else:
    namelist=dat["System Number"][istart:]

for i,sysi in enumerate(namelist):
    f=open("sb9_position.txt","a")

    ##### GET positions #####
#    if True:
    try:
        pradec=dat["2000.0 coordinates"][sysi]
        c=SkyCoord.from_name("J"+pradec, parse=True)
        ra=c.ra.degree
        dec=c.dec.degree

#    if True:
        width = u.Quantity(5, u.arcsec)
        height = u.Quantity(5, u.arcsec)
        #GAIA
        r = Gaia.query_object_async(coordinate=c, width=width, height=height)
        plx=None
        if len(r["parallax"]) == 0:
            sw = False
        elif type(r["parallax"][0]) == np.float64:
            plx=r["parallax"][0]
            sw = True
        else:
            sw = False
        print("GAIA",plx)

        Simbad.SIMBAD_URL = "http://simbad.u-strasbg.fr/simbad/sim-script"
        Simbad.add_votable_fields("parallax","flux(V)","flux(R)","flux(J)","flux(H)","flux(K)")

        result_table = Simbad.query_region(c, radius='0d0m5s')
        print(result_table)
        if result_table is None:
            plxs=np.nan
            magcom="|||||"            
        elif len(result_table) == 1:
            plxs=result_table["PLX_VALUE"].item()
            V=result_table["FLUX_V"].item()
            R=result_table["FLUX_R"].item()        
            J=result_table["FLUX_J"].item()
            H=result_table["FLUX_H"].item()
            K=result_table["FLUX_K"].item()            
            magcom="|"+str(V)+"|"+str(R)+"|"+str(J)+"|"+str(H)+"|"+str(K)
        else:
            plxs=result_table["PLX_VALUE"][0]
            V=result_table["FLUX_V"][0]
            R=result_table["FLUX_R"][0] 
            J=result_table["FLUX_J"][0]
            H=result_table["FLUX_H"][0]
            K=result_table["FLUX_K"][0]     
            magcom="|"+str(V)+"|"+str(R)+"|"+str(J)+"|"+str(H)+"|"+str(K)

        #eplx=result_table["PLX_ERROR"].item()
        if plxs == plxs:
            f.write(str(sysi)+"|"+str(ra)+"|"+str(dec)+"|"+str(plxs)+"|"+str(plx)+magcom+"\n")
        else:
            f.write(str(sysi)+"|"+str(ra)+"|"+str(dec)+"|None|"+str(plx)+magcom+"\n")

    except:
        try:
            pradec=dat["2000.0 coordinates"][sysi]
            c=SkyCoord.from_name("J"+pradec, parse=True)
            ra=c.ra.degree
            dec=c.dec.degreeV
            f.write(str(sysi)+"|"+str(ra)+"|"+str(dec)+"|||||||"+"\n")
        except:
            f.write(str(sysi)+"|||||||||"+"\n")

    f.close()




#dat3=pd.read_csv("../database/sb9/Alias.dta",delimiter="|",dtype={"System number":"int","Catalog name":"str","ID in that catalog":"str"})

#dat.to_pickle("../database/pickles/sb9.Main.pickle")
#dat2.to_pickle("../database/pickles/sb9.Orbits.pickle")
#dat3.to_pickle("../database/pickles/sb9.Alias.pickle")

