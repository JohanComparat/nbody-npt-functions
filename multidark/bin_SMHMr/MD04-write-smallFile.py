from MultiDark import * 
import os
import numpy as n
import glob

snList = n.array(glob.glob(os.path.join(os.environ["MD04"], "snapshots", "hlist_*.list")))
print snList
box = MultiDarkSimulation(Lbox=400.0 * uu.Mpc, boxDir = "MD_0.4Gpc", snl =snList   ,zsl = None, zArray = n.arange(0.2,2.4,1e-1), Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))
box.columnDict = box.columnDictHlist

t0=time.time()

for ii in n.arange(len(box.snl)):
	box.writeSAMcatalog(ii, mmin=100*box.Melement)

print time.time()-t0