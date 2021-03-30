import numpy as np
import sys
def gen(n,i,a,b):
    f=open("data/test"+str(i+5),"w")
    randnums= np.random.randint(a,b,n)
    for e in randnums:
        f.write(str(e)+" ")
    f.write("\n")
    f.close()
for n in range(2,int(sys.argv[1])):
    gen(2**n,n,-n*n,n*n)