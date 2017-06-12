from MultiDark import *
box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_2.5Gpc",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" , "hlist_?.?????.list"))) ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)


for ii in n.arange(len(box.snl)):
	box.computeDoubleDistributionFunction(ii,'Vpeak', 'mvir', n.arange(0,3.5,0.01), n.arange(8,16,0.01))
	box.combinesDoubleDistributionFunction(ii,'Vpeak', 'mvir', n.arange(0,3.5,0.01), n.arange(8,16,0.01),type = "Central")
	box.combinesDoubleDistributionFunction(ii,'Vpeak', 'mvir', n.arange(0,3.5,0.01), n.arange(8,16,0.01),type = "Satellite")

