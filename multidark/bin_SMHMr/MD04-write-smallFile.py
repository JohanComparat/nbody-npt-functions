import MultiDark as md
import os
import numpy as n
import glob

snList = n.array(glob.glob(os.path.join(os.environ["MD04"], "snapshots", "out_*.list")))

box = md.MultiDarkSimulation(Lbox=400.0 * uu.Mpc, boxDir = "MD_0.4Gpc", snl =snList   ,zsl = None, zArray = n.arange(0.2,2.4,1e-1), Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writeSAMcatalog(ii, mmin=1000*box.Melement)

