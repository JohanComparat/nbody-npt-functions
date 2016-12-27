from MultiDark import *

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc_new_rockS")

all=n.array([ 0.50320, 0.403200, 0.600800]) #, 0.89510, 0.74980 ])
for ii in range(len(all)):
	a= all[ii]
	ll = n.array( glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/*"+str(a)+"*.fits" ) )
	#box.compute2PCF(ll, rmax=140, Nmax=4000000, vmin=65, dr = 2., name="rmax_140")
	#box.compute2PCF(ll, rmax=15, Nmax=4000000, vmin=65, dr = 0.1, name="rmax_015")
	box.compute2PCF(ll, rmax=50, Nmax=4000000, vmin=65, dr = 1., name="rmax_050")

