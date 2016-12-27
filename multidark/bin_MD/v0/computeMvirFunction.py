from MultiDark import *
box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_2.5Gpc",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots" , "hlist_?.?????.list"))) ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)

bins = n.arange(8,16,0.01)
quantity = 'mvir'
for ii in range(len(box.snl)):
	print box.snl[ii]
	box.computeSingleDistributionFunction(ii, quantity, bins)

output_dir = join(box.wdir,box.boxDir,"properties",quantity)
n.savetxt(join(output_dir , quantity+ '.bins'), n.transpose([bins]))