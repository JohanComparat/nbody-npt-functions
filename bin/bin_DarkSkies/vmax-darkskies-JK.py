from DarkSkies import *

box = DarkSkiesSimulation()

for ii in n.arange(len(box.snl)):
	print box.snl[ii]
	fileList = glob.glob(box.snl[ii][:-4]+"*_PM_Nb_*.fits")
	rootname = os.path.basename(box.snl[ii])[:-4]
	box.computeSingleDistributionFunctionJKresampling( fileList, rootname, "vmax", 10**n.arange(0,3.5,0.01),  Ljk = 400., overlap = 1. )

