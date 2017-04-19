from MultiDark import *
import sys
import glob
import os 



snList= n.array(glob.glob(os.path.join(os.environ["MD10"], "snapshots", "out_*.list")))
print snList
box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

t0=time.time()

for ii in n.arange(len(box.snl)):
	box.writeSAMcatalog(ii, mmin=100*box.Melement)

print time.time()-t0

sys.exit()

snList= n.array(glob.glob(os.path.join(os.environ["MD10"], "snapshots", "hlist_*.list")))
print snList
box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))
box.columnDict = box.columnDictHlist

t0=time.time()

for ii in n.arange(len(box.snl)):
	box.writeSAMcatalog(ii, mmin=100*box.Melement)

print time.time()-t0