from MultiDark import *

snList= n.array([ "/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_1.00000.list" , 
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.89510.list",
"/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/snapshots/hlist_0.74980.list"])

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc_new_rockS",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

zs= n.array([0, 0.11719360965255277, 0.33368898372899425])
bd = (zs[1:] + zs[:-1])/2.
ds = c2.comoving_distance(bd) * c2.h
ts = c2.lookback_time(zs)
td = c2.lookback_time(bd)
dmins = n.hstack(([0.], ds ))
dmaxs = n.hstack((ds, [ 1000.]))

ii=1
box.cornerLCpositionCatalog(ii, DMIN=dmins[ii], DMAX=dmaxs[ii], vmin = 190)

