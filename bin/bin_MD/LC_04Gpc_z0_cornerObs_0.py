from MultiDark import *

snList=  n.array(['/data2/DATA/eBOSS/Multidark-lightcones/MD_0.4Gpc/snapshots/hlist_1.00000.list','/data2/DATA/eBOSS/Multidark-lightcones/MD_0.4Gpc/snapshots/hlist_0.90430.list',])

box = MultiDarkSimulation(Lbox=400.0 , boxDir = "MD_0.4Gpc",snl =snList )


zs= n.array([0, 0.10582771204246377 ])
bd = (zs[1:] + zs[:-1])/2.
ds = c2.comoving_distance(bd) * c2.h
ts = c2.lookback_time(zs)
td = c2.lookback_time(bd)
dmins = n.hstack(([0.], ds ))
dmaxs = n.hstack((ds, [ 400.]))

ii=0
box.cornerLCpositionCatalog(ii, DMIN=dmins[ii], DMAX=dmaxs[ii], vmin = 65)
