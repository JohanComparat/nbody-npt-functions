from MultiDark import *

box = MultiDarkSimulation(Lbox=400.0 * uu.Mpc, boxDir = "MD_0.4Gpc")

all=n.array([ 1.00000, 0.50000, 0.31800, 0.40700, 0.59240, 0.74230, 0.90430 ])
for ii in range(len(all)):
	a= all[ii]
	ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_0.4Gpc/snapshots/*"+str(a)+"*.fits" ) )
	box.compute2PCF(ll, rmax=140, Nmax=4000000, vmin=65, dr = 2., name="rmax_140")
	box.compute2PCF(ll, rmax=15, Nmax=4000000, vmin=65, dr = 0.1, name="rmax_015")

