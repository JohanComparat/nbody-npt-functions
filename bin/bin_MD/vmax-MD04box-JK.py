from MultiDark import *

snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_0.4Gpc", "snapshots", "out_*.list")))
box = MultiDarkSimulation(Lbox=400.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_0.4Gpc", snl = snl ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 96300000.0)

for ii in n.arange(len(box.snl)):
	fileList = glob.glob(box.snl[ii][:-5]+"*_PM_Nb_*.fits")
	rootname = box.snl[ii][:-5].split('/')[-1]
	box.computeSingleDistributionFunctionJKresampling( fileList, rootname, "vmax", 10**n.arange(0,3.5,0.01),  Ljk = 40., overlap = 1. )

