import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import sys


dat4=pd.read_csv("../siwiyn/siwiyn_position_clean2.txt",delimiter="|",dtype={"System Number":"int","Name":"str","ra":"str","dec":"str","Simbad plx":"float","GAIA plx":"float","V":"float","R":"float","J":"float","H":"float","K":"float"})
dat4=dat4.set_index("System Number")

datx=pd.read_csv("../siwiyn/siwiyn.tsv",delimiter="|",comment="#")
mask=~datx.duplicated(subset='Name')

#remove dupicated sysnum
newd=datx[mask]
#newd["radec"]= ""
newd["Simbad plx"]= np.nan
newd["GAIA plx"]= np.nan
newd["V"]= np.nan
newd["R"]= np.nan
newd["J"]= np.nan
newd["H"]= np.nan
newd["K"]= np.nan
#newd=newd.set_index("Name")

#print(newd)
#sys.exit()

for i,ind in enumerate(newd.index):
    print(i,ind)
    name=newd["Name"][ind]
    if True:
        
#    try:
        maskx=dat4["Name"]==name
        for comp in ["V","R","J","H","K","Simbad plx","GAIA plx"]:
            if len(dat4[comp].values[maskx])>0:
                newd[comp][ind]=dat4[comp].values[maskx][0]

#            newd[]
#        newd.loc[sysi,"name"]=dat["Name"].values[0]
        
        #print(i,sysi)
#        newd.loc[sysi,"Sp1"]=dat["Spectral type component 1"][i]#.values[maskx][0]
#        newd.loc[sysi,"Sp2"]=dat["Spectral type component 2"][i]#.values[maskx][0]        
#        newd.loc[sysi,"radec"]=dat["2000.0 coordinates"][i]#.values[maskx][0]
#        newd.loc[sysi,"Simbad plx"]=dat4.loc[sysi,"Simbad plx"]
#        newd.loc[sysi,"GAIA plx"]=dat4.loc[sysi,"GAIA plx"]
#        newd.loc[sysi,"V"]=dat4.loc[sysi,"V"]
#        newd.loc[sysi,"R"]=dat4.loc[sysi,"R"]
#        newd.loc[sysi,"J"]=dat4.loc[sysi,"J"]
#        newd.loc[sysi,"H"]=dat4.loc[sysi,"H"]
#        newd.loc[sysi,"K"]=dat4.loc[sysi,"K"]
#    except:
#        print("No data for ",sysi)
    

newd.to_csv("siwiyn.csv")
