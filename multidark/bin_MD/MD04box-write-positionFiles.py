from MultiDark import *

snList = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_0.4Gpc/snapshots/out_*.list"))

box = MultiDarkSimulation(Lbox=400.0 * uu.Mpc, boxDir = "MD_0.4Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))


for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin=30, mmin=1*10**9)


