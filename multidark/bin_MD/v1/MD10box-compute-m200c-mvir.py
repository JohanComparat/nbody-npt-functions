from MultiDark import *
box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_1Gpc_new_rockS",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_1Gpc_new_rockS", "snapshots", "hlist_*.list"))) ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 1510000000.0)


for ii in n.arange(len(box.snl)):
	box.computeDoubleDistributionFunction(ii,'M200c', 'mvir',n.arange(8,16,0.01), n.arange(8,16,0.01))
	box.combinesDoubleDistributionFunction(ii,'M200c', 'mvir',n.arange(8,16,0.01), n.arange(8,16,0.01),type = "Central")
	box.combinesDoubleDistributionFunction(ii,'M200c', 'mvir',n.arange(8,16,0.01), n.arange(8,16,0.01),type = "Satellite")

