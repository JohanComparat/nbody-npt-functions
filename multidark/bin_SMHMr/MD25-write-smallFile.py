from MultiDark import *
import sys
import glob
import os
import time


snList= n.array(glob.glob(os.path.join(os.environ["MD25"], "snapshots", "out_*.list")))
print snList 
box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_2.5Gpc",snl =  snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)
#box.columnDict = box.columnDictHlist25

t0=time.time()

for ii in n.arange(len(box.snl)):
	box.writeSAMcatalog(ii, mmin=100*box.Melement)

print time.time()-t0

sys.exit()

snList= n.array(glob.glob(os.path.join(os.environ["MD25"], "snapshots", "hlist_*.list")))
print snList 
box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_2.5Gpc",snl =  snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)
box.columnDict = box.columnDictHlist25

t0=time.time()

for ii in n.arange(len(box.snl)):
	box.writeSAMcatalog(ii, mmin=100*box.Melement)

print time.time()-t0

