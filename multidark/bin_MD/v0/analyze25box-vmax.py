from MultiDark import *
box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_2.5Gpc",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" , "hlist_?.?????.list"))) ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)

#bins = n.loadtxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/vmax/vmax.bins",unpack=True)

for ii in n.arange(len(box.snl)):
    box.computeSingleDistributionFunction(ii,'vmax', n.arange(0,3.5,0.01))
    box.combinesSingleDistributionFunction(ii,'vmax', n.arange(0,3.5,0.01),type = "Central")
    box.combinesSingleDistributionFunction(ii,'vmax', n.arange(0,3.5,0.01),type = "Satellite")

import sys
sys.exit()


names=n.array(["M200b","M200c","M2500c","M500c","Macc","mvir"])
for qty in names :
    for ii in n.arange(len(box.snl)):
        box.combinesSingleDistributionFunction(ii, qty, 10**n.arange(8,16,0.01),type = "Central")
        box.combinesSingleDistributionFunction(ii, qty, 10**n.arange(8,16,0.01),type = "Satellite")

names=n.array(["Vacc","vmax"])
for qty in names :
    for ii in n.arange(len(box.snl)):
        box.combinesSingleDistributionFunction(ii, qty, 10**n.arange(0,4.5,0.01),type = "Central")
        box.combinesSingleDistributionFunction(ii, qty, 10**n.arange(0,4.5,0.01),type = "Satellite")
