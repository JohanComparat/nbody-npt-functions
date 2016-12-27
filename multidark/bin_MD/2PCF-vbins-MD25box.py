from MultiDark import *

box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc, boxDir = "MD_2.5Gpc")

all=n.array([   46, 11, 10, 7,9 ]) # 80,  22, 74,
for ii in range(len(all)):
	a= all[ii]
	ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_"+str(a)+"_Nb*.fits" ) )
	box.compute2PCF(ll, rmax=50, Nmax=4000000, vmin=65, dr = 1., name="rmax_050")
	#box.compute2PCF(ll, rmax=140, Nmax=4000000, vmin=65, dr = 2., name="rmax_140")
	#box.compute2PCF(ll, rmax=15, Nmax=4000000, vmin=65, dr = 0.1, name="rmax_015")

