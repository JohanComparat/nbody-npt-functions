from DarkSkies import *

box = DarkSkiesSimulation()


for ii in n.arange(len(box.snl)):
	print box.snl[ii]
	fileList = glob.glob(box.snl[ii][:-4]+"*_PM_Nb_*.fits")
	rootname = os.path.basename(box.snl[ii])[:-4]
	box.computeSingleDistributionFunctionJKresampling( fileList, rootname, "mvir", n.arange(8,16,0.05),  Ljk = 400., overlap = 1. )

