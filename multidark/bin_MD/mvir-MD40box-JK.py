from MultiDark import *
box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_4Gpc",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_4Gpc", "snapshots", "out_128.list"))) ,zsl = None,zArray = n.arange(0.,4.,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 96000000000.0)
box.snl.sort()

for ii in n.arange(len(box.snl)):
	print box.snl[ii]
	fileList = glob.glob(box.snl[ii][:-5]+"*_PM_Nb_*.fits")
	rootname = box.snl[ii][:-5].split('/')[-1]
	box.computeSingleDistributionFunctionJKresampling( fileList, rootname, "mvir", n.arange(8,16,0.05),  Ljk = 400., overlap = 1. )

