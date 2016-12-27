from MultiDark import *
import glob 
snList=  n.array(glob.glob("/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_*.list"))

snList=  n.array(["/home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_4GpcNW/snapshots/out_16.list"])

box = MultiDarkSimulation(Lbox=4000.0 * uu.Mpc, boxDir = "MD_4GpcNW",snl = snList ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 200., mmin=9*10**11)


