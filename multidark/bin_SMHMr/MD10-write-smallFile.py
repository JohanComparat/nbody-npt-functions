from MultiDark import *
import sys
import glob
import os 

snList= n.array(glob.glob(os.path.join(os.environ["MD10"], "snapshots", "out_*.list")))

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writeSAMcatalog(ii, mmin=1000*box.Melement)

