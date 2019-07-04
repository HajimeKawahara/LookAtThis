import numpy as np
g = open("sb9_position_clean.txt","w")

f = open("sb9_position.txt","r")
for line in f:
    arr=(line.split("|"))
    #['3090', '', '', '', '', '', '', '', '', '\n']
    l=len(arr[1])+len(arr[2])+len(arr[3])+len(arr[4])+len(arr[5])+len(arr[6])+len(arr[5])+len(arr[6])
    if l == 0:
        print(line)
    else:
        g.write(line)
f.close()
g.close()
