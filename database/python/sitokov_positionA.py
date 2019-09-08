import pandas as pd
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import numpy as np
import argparse
import sys
from time import sleep

parser = argparse.ArgumentParser(description='TOKOVININ19')
parser.add_argument('-i', nargs=1, default=[0], help='start index', type=int)
parser.add_argument('-f', nargs=1, help='file', type=str)

args = parser.parse_args()    
istart = args.i[0]

dat=pd.read_csv("../sitok/tokovinin19_A.txt",delimiter="&",comment="#")
#GAIA distance
#f=open("siwiyn_position.txt","a")
#f.write("System Number|name|RA(2000)|DEC(2000)|Simbad plx|GAIA plx|V|R|J|H|K"+"\n")
#f.close()
if args.f:
    namelist=np.loadtxt(args.f[0],dtype=int)
else:
    namelist=range(istart,len(dat))

for i,sysi in enumerate(namelist):
    gw=open("namesimbad.txt","a")
    f=open("sitok_positionA.txt","a")
    name=dat["Name"].values[sysi]
    print(i,name)
    
    sleep(1.5)
    if True:
#    try:
        radec=dat["pos"][sysi]
        sep=dat["sep"][sysi]
        sp=dat["sp"][sysi]

        print(radec)
        c = SkyCoord("J"+radec, unit=(u.hourangle, u.deg))
        
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
        Simbad.add_votable_fields("parallax","flux(V)","flux(R)","flux(J)","flux(H)","flux(K)","pmra","pmdec")

        result_table = Simbad.query_region(c, radius='0d0m5s')
        print(result_table)
        if result_table is None:
            plxs=np.nan
            magcom="|||||"            
        else:
            simbadname=result_table["MAIN_ID"][0]
            print("--------------------------")
            print(name,simbadname)
            gw.write(name+"|"+simbadname.decode("utf-8")+"\n")
            print("--------------------------")
            name=simbadname.decode("utf-8") #rename
            
            plxs=result_table["PLX_VALUE"][0]
            V=result_table["FLUX_V"][0]
            R=result_table["FLUX_R"][0] 
            J=result_table["FLUX_J"][0]
            H=result_table["FLUX_H"][0]
            K=result_table["FLUX_K"][0]     
            magcom="|"+str(V)+"|"+str(R)+"|"+str(J)+"|"+str(H)+"|"+str(K)
            rasimbad=result_table["RA"][0]
            decsimbad=result_table["DEC"][0]
            pmra=result_table["PMRA"][0]
            pmdec=result_table["PMDEC"][0]
            radecinfo="|"+str(rasimbad)+"|"+str(decsimbad)+"|"+str(pmra)+"|"+str(pmdec)
        #eplx=result_table["PLX_ERROR"].item()
        if plxs == plxs:
            f.write(str(sysi)+"|"+name+"|"+str(radec)+"|"+str(sep)+"|"+str(sp)+"|"+str(plxs)+"|"+str(plx)+magcom+"|"+radecinfo+"\n")
        else:
            f.write(str(sysi)+"|"+name+"|"+str(radec)+"|"+str(sep)+"|"+str(sp)+"|None|"+str(plx)+magcom+"|"+radecinfo+"\n")

#    except:
#        try:
#            radec=dat["pos"][sysi]
#            sep=dat["sep"][sysi]
#            c = SkyCoord(radec, unit=(u.hourangle, u.deg))
#            f.write(str(sysi)+"|"+str(radec)+"|"+str(sep)+"|"+str(sp)+"||||||"+"\n"#)
#        except:
#            f.write(str(sysi)+"||||||||||"+"\n")

    f.close()
    gw.close()




