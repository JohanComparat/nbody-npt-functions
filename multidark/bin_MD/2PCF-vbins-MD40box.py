from MultiDark import *

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4Gpc")

all=n.array([ 128, 97, 105, 66, 79, 82, 112, 124, 87, 91, 97 ])
for ii in range(len(all)):
	a= all[ii]
	ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/hlist_"+str(a)+"_Nb*.fits" ) )
	#box.compute2PCF(ll, rmax=140, Nmax=4000000, vmin=65, dr = 2., name="rmax_140")
	box.compute2PCF(ll, rmax=50, Nmax=4000000, vmin=65, dr = 1., name="rmax_050")
	#box.compute2PCF(ll, rmax=15, Nmax=4000000, vmin=65, dr = 0.1, name="rmax_015")


all=n.array([ 97, 112, 124, 87, 91, 97 ])
for ii in range(len(all)):
	a= all[ii]
	ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/hlist_"+str(a)+"_Nb*.fits" ) )
	box.compute2PCF(ll, rmax=15, Nmax=4000000, vmin=65, dr = 0.1, name="rmax_015")
	box.compute2PCF(ll, rmax=140, Nmax=4000000, vmin=65, dr = 2., name="rmax_140")
	#box.compute2PCF(ll, rmax=50, Nmax=4000000, vmin=65, dr = 1., name="rmax_050")

