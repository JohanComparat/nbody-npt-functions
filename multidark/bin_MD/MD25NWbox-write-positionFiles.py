from MultiDark import *
snList =  n.array(["/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_2.5GpcNW/snapshots/out_80.list"])

box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc, boxDir = "MD_2.5GpcNW",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 100., mmin=2*10**11)

