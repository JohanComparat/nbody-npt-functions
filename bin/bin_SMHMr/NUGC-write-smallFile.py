from MultiDark import *
import sys
import glob
import os
import astropy.io.fits as fits

summary_file_name = os.path.join(os.environ["NUGC_DIR"], 'output_NUGC.fits')

sf = fits.open(summary_file_name)[1].data

t0=time.time()


snap = sf[2]
print snap

snList= n.array([os.path.join(os.environ["NUGC_DIR"], "snapshots", snap['snap_name'])])

box = MultiDarkSimulation(Lbox=snap['box_length'] * uu.Mpc, boxDir = "NUGC", snl = snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 68.0 * uu.km / (uu.s * uu.Mpc), Melement= snap['mass_particle'])


if snap['N_columns']==39 :
	box.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'c_to_a_500c': 33, 'Ax_500c': 34, 'Ay_500c': 35, 'Az_500c': 36, 'TU': 37, 'pid': 38}
	

if snap['N_columns']==43 :
	box.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'c_to_a_500c': 33, 'Ax_500c': 34, 'Ay_500c': 35, 'Az_500c': 36, 'TU': 37, 'M_pe_Behroozi': 38, 'M_pe_Diemer': 39, 'Halfmass_Radius': 40, 'rvmax': 41, 'pid': 42}
	
box.writeSAMcatalog(0, mmin=300*box.Melement)

print time.time()-t0

