import pandas as pd
import numpy as np
import typetomass as ttm

if __name__ == "__main__":
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

        #print(sp1,mass1,sp2,mass2)
        #print("Mtot=",totalmass,"P=",P,"(d) a=",a,"au")
        #print("d=",1/(plx_simbad/1000),"pc")
        #print("theta gaia=",theta_gaia,"mas","theta simbad=",theta_simbad,"mas")
        if theta_simbad > 100:
            print(sp1,mass1,sp2,mass2,"theta simbad=",theta_simbad,"mas")
#            print("Mtot=",totalmass,"P=",P,"(d) a=",a,"au")
#            print("d=",1/(plx_simbad/1000),"pc")
#            print("theta gaia=",theta_gaia,"mas","theta simbad=",theta_simbad,"mas")
            
