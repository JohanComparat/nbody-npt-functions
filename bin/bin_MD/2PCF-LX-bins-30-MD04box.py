from MultiDark import *

box = MultiDarkSimulation(Lbox=400.0 * uu.Mpc, boxDir = "MD_0.4Gpc")

ll = n.array( glob.glob( "/home/comparat/darksim/MD/MD_0.4Gpc/catalogs//out_54*.fits" ) )
box.compute2PCF_LX(ll, rmax=30, dr = 0.1, Nmax=1000000, vmin=41.5, dlogBin=0.2,  name="rmax_30")

