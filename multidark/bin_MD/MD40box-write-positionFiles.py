from MultiDark import *

snList =  n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/out_12?.list"))

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4Gpc",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 200., mmin=9*10**11)


snList =  n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/out_5?.list"))

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4Gpc",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 200., mmin=9*10**11)


snList =  n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/out_6?.list"))

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4Gpc",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 200., mmin=9*10**11)


snList =  n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/out_7?.list"))

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4Gpc",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 200., mmin=9*10**11)


snList =  n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/out_9?.list"))

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4Gpc",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 200., mmin=9*10**11)

