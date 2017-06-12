from MultiDark import *
import sys
import glob
import os
import astropy.io.fits as fits

sf = fits.open(os.path.join(os.environ["MD10"], "output_MD_1.0Gpc.fits"))[1].data

snList= n.array([os.path.join(os.environ["MD10"], "snapshots", "out_"+name+".list") for name in sf['snap_name']])

box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

t0=time.time()

for ii in n.arange(len(box.snl)):
	print ii, box.snl[ii], os.path.isfile(box.snl[ii][:-5]+"_cluster_Nb_0.fits")
	if os.path.isfile(box.snl[ii][:-5]+"_cluster_Nb_0.fits")==False :
		if sf['N_columns'][ii]==34 :
			box.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'pid': 32}
		if sf['N_columns'][ii]==41 :
			box.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'c_to_a_500c': 33, 'Ax_500c': 34, 'Ay_500c': 35, 'Az_500c': 36, 'TU': 37, 'M_pe_Behroozi': 38, 'M_pe_Diemer': 39, 'pid': 40}
			
		box.writeSAMcatalog(ii, mmin=300*box.Melement)

print time.time()-t0


from MultiDark import *
import sys
import glob
import os 



snList= n.array(glob.glob(os.path.join(os.environ["MD10"], "snapshots", "out_*.list")))
print snList
box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))

t0=time.time()

for ii in n.arange(len(box.snl)):
	box.writeSAMcatalog(ii, mmin=100*box.Melement)

print time.time()-t0

sys.exit()

snList= n.array(glob.glob(os.path.join(os.environ["MD10"], "snapshots", "hlist_*.list")))
print snList
box = MultiDarkSimulation(Lbox=1000.0 * uu.Mpc, boxDir = "MD_1Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc))
box.columnDict = box.columnDictHlist

t0=time.time()

for ii in n.arange(len(box.snl)):
	box.writeSAMcatalog(ii, mmin=100*box.Melement)

print time.time()-t0