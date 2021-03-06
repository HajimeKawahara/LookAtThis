import pandas as pd
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import numpy as np
import argparse
import sys
from time import sleep

parser = argparse.ArgumentParser(description='SIXTH: get information from Simbad and GAIA TAP (DR1)')
parser.add_argument('-i', nargs=1, default=[0], help='start index', type=int)
parser.add_argument('-f', nargs=1, help='file', type=str)

args = parser.parse_args()    
istart = args.i[0]

datx=pd.read_csv("../siwiyn/siwiyn.tsv",delimiter="|",comment="#")
#_RAJ2000|_DEJ2000|WDS|Name|HD|HIP|Date|PA|Sep|Dmag|Wave|FWHM|f_FWHM|SimbadName|_RA|_DE
mask=~datx.duplicated(subset='Name')
dat=datx[mask]

#dat=dat.set_index("Name")

#print(dat["Name"].values[0])
print(len(datx),len(dat))
#sys.exit()

#GAIA distance
#f=open("siwiyn_position.txt","a")
#f.write("System Number|name|RA(2000)|DEC(2000)|Simbad plx|GAIA plx|V|R|J|H|K"+"\n")
#f.close()
if args.f:
    namelist=np.loadtxt(args.f[0],dtype=int)
else:
    namelist=range(istart,len(dat))

for i,sysi in enumerate(namelist):
    f=open("siwiyn_position.txt","a")
    name=dat["Name"].values[sysi]
    print(i,name)
    
    sleep(1)
#    if True:
    try:
        ra=dat["_RAJ2000"][sysi]
        dec=dat["_DEJ2000"][sysi]
        c = SkyCoord(ra+" "+dec, unit=(u.hourangle, u.deg))
        
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
            f.write(str(sysi)+"|"+name+"|"+str(ra)+"|"+str(dec)+"|"+str(plxs)+"|"+str(plx)+magcom+"\n")
        else:
            f.write(str(sysi)+"|"+name+"|"+str(ra)+"|"+str(dec)+"|None|"+str(plx)+magcom+"\n")

    except:
        try:
            ra=dat["_RAJ2000"][sysi]
            dec=dat["_DEJ2000"][sysi]
            c = SkyCoord(ra+" "+dec, unit=(u.hourangle, u.deg))
            f.write(str(sysi)+"|"+str(ra)+"|"+str(dec)+"|||||||"+"\n")
        except:
            f.write(str(sysi)+"|||||||||"+"\n")

    f.close()





