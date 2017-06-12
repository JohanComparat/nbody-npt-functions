from MultiDark import *
import sys
import glob

snList= n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc/snapshots/out_*.list"))

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii,  vmin=50, mmin=2*10**10)

