import numpy as np
g = open("siwiyn_position_clean2.txt","w")

f = open("siwiyn_position_clean.txt","r")
for line in f:
    arr=(line.split("|"))
    #['3090', '', '', '', '', '', '', '', '', '\n']
    l=len(arr[3])+len(arr[4])+len(arr[5])+len(arr[6])+len(arr[5])+len(arr[6])
    if l == 0:
        print(arr[0])
    else:
        g.write(line)
f.close()
g.close()
