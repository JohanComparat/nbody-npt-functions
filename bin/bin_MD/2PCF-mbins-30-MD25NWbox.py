from MultiDark import *

box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc, boxDir = "MD_2.5GpcNW")

ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5GpcNW/snapshots/out_80*.fits" ) )
box.compute2PCF_MASS(ll, rmax=30, dr = 0.1, Nmax=1000000, vmin=n.log10(box.Melement)+2, dlogBin=0.05,  name="rmax_30")


