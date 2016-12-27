from MultiDark import *
box = MultiDarkSimulation(Lbox=400.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_0.4Gpc",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_0.4Gpc", "snapshots", "hlist_?.?????.list"))) ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 96300000.0)

for ii in n.arange(len(box.snl)):
	box.computeDoubleDistributionFunction(ii,'Vpeak', 'mvir', n.arange(0,3.5,0.01), n.arange(8,16,0.01))
	box.combinesDoubleDistributionFunction(ii,'Vpeak', 'mvir', n.arange(0,3.5,0.01), n.arange(8,16,0.01),type = "Central")
	box.combinesDoubleDistributionFunction(ii,'Vpeak', 'mvir', n.arange(0,3.5,0.01), n.arange(8,16,0.01),type = "Satellite")


import sys
sys.exit()

nameA, nameB, binsA, binsB = 'Vpeak', 'mvir', n.arange(0,3.5,0.01), n.arange(8,16,0.01)
type = "Central" 
output_dir = join(box.wdir,box.boxDir,"properties",nameA+"-"+nameB)
nameSnapshot = box.snl[ii].split('/')[-1][:-5]
pklList = n.array(glob.glob(join(output_dir, nameSnapshot + "_" + nameA+"-"+nameB +"_"+type+"_*.pkl")))

nnA = n.empty( [len(pklList),len(binsA)-1] ) 
nnB = n.empty( [len(pklList),len(binsB)-1] ) 
dataAB = n.empty( [len(pklList),len(binsA)-1,len(binsB)-1] ) 
for jj in range(len(pklList)):
    f=open(pklList[jj],'r')
    nnAinter, nnBinter, dataABinter = cPickle.load(f)
    nnA[jj] = nnAinter
    nnB[jj] = nnBinter
    dataAB[jj] = dataABinter[0]
    f.close()

n.savetxt(join(output_dir,"hist-"+type+"-"+nameA+"-"+nameSnapshot[6:]+".dat"),n.transpose([binsA[:-1], binsA[1:], nnA.sum(axis=0)]))
n.savetxt(join(output_dir,"hist-"+type+"-"+nameB+"-"+nameSnapshot[6:]+".dat"),n.transpose([binsB[:-1], binsB[1:], nnB.sum(axis=0)]))
n.savetxt(join(output_dir, "hist2d-"+type+"-"+ nameA+"-"+nameB + "-"+ nameSnapshot[6:] + ".dat"), dataAB.sum(axis=0))

