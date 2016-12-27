from MultiDark import *

snl = n.array([ 
join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_10.list"), 
join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_19.list"), 
join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_74.list"), 
join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_80.list"), 
join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_11.list"), 
join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_22.list"), 
join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_46.list"), 
join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_7.list"), 
join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_9.list") ])

box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_2.5Gpc",snl =  snl ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)


for ii in n.arange(len(box.snl)):
	box.writePositionCatalogPM(ii, vmin = 100., mmin=2*10**11)


snl = n.array([ join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc", "snapshots" ,"out_30.list") ])

box = MultiDarkSimulation(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_2.5Gpc",snl =  snl ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)

box.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'c_to_a_500c': 33, 'Ax_500c': 34, 'Ay_500c': 35, 'Az_500c': 36, 'TU': 37, 'M_pe_Behroozi': 38, 'M_pe_Diemer': 39, 'Halfmass_Radius': 40, 'pid': 41}

for ii in n.arange(len(box.snl)):
	fileList = glob.glob(box.snl[ii][:-5]+"*_PM_Nb_*.fits")
	rootname = box.snl[ii][:-5].split('/')[-1]
	box.computeSingleDistributionFunctionJKresampling( fileList, rootname, "mvir", n.arange(8,16,0.1),  Ljk = 250., overlap = 1. )
