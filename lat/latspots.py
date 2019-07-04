import pandas as pd
import numpy as np
import typetomass as ttm

#Seq|Star|RAJ2000|DEJ2000|Dist|n_Dist|Age|Hmag|SpType|MA|l_MB|MB|l_rho|rho|l_Per|Per|x_Per|l_Ecc|Ecc|l_acrit|acrit|Ref|SimbadName

def is_float_expression(s):
    try:
        val = float(s)
        return True
    except ValueError:
        return False
    
if __name__ == "__main__":
    dat=pd.read_csv("../database/spots/spots.tsv",delimiter="|")
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
        if theta > 100.0:
            print(name,sp,mass1,mass2,"theta=",theta)
