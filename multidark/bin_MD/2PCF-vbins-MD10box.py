from MultiDark import *

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc")

ll = n.array( glob.glob( "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc/snapshots/out_128*.fits" ) )
box.compute2PCF_MASS(ll, rmax=140, dr = 0.1, Nmax=2000000, vmin=n.log10(box.Melement)+1, dlogBin=0.1,  name="rmax_140")

