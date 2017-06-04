from MultiDark import *
import sys
import glob
import os

snList= n.array(glob.glob(os.path.join(os.environ["MD10"], "snapshots", "out_*p.list")))
snList.sort()

print snList
box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

t0=time.time()

for ii in n.arange(4,len(box.snl),1):
        print ii, box.snl[ii]
        box.writeCLUSTERcatalog(ii, mmin=6*10**12)

print time.time()-t0
