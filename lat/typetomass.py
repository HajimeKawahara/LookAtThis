
#very rough estimater of mass
def type2mass(Sp,dbV,unknownmass=0.0):
    if Sp=="sdB":
        mass=0.6 # white dwarf
        return mass

    #main squence
    try:
        MS=None
        if len(Sp) == 2:
            MS=Sp[0:2]
        elif Sp[2] != "I":
            MS=Sp[0:2]
        if MS is not None:
            mask=(dbV["Sp"]==MS)
            if len(dbV["mass"][mask]) > 0:
                mass=dbV["mass"][mask].values[0]
                return mass
            
        return unknownmass

    except:
        return unknownmass

    
if __name__ == "__main__":
    import pandas as pd
    from astroquery.simbad import Simbad
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia
    import numpy as np
    
    #    dat=pd.read_csv("../database/sb9/Main.dta",delimiter="|",dtype={"System Number":"int","1900.0 coordinates":"str","2000.0 coordinates":"str","Component":"str","Magnitude of component 1":"float","Filter component 1":"str","Magnitude of component 2":"float","Filter component 2":"str","Spectral type component 1":"str","Spectral type component 2":"str"})
    
    dat=pd.read_csv("../database/sb9/sb9.csv")

    dbV=pd.read_csv("../database/sp/spV.txt",delimiter=",")
#    print(dbV["mass"][dbV["Sp"]=="O8"].values[0])
    for sp in dat["Sp1"]:
        mass1=type2mass(sp,dbV,unknownmass=1.0)
        print(sp,mass)
