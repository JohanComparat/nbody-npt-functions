from MultiDark import *
import sys
import glob
import os
import astropy.io.fits as fits

ii = int(sys.argv[1])
print ii

sf = fits.open(os.path.join(os.environ["MD10"], "output_MD_1.0Gpc.fits"))[1].data
snList= n.array([os.path.join(os.environ["MD10"], "snapshots", "out_"+name+".list") for name in sf['snap_name']])
#snList= n.array([os.path.join(os.environ["MD10"], "snapshots", "out_99p.list")])
box = MultiDarkSimulation(Lbox=1000.0 * u.Mpc, boxDir = "MD_1Gpc",snl =snList   ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * u.km / (u.s * u.Mpc))
t0=time.time()

#sf['snap_name'][ii]

name = os.path.join(os.environ["MD10"], "snapshots", "out_"+sf['snap_name'][ii]+".list")

#for ii in n.arange(len(box.snl)):
print sf[ii]
if sf['N_columns'][ii]==34 :
	box.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'pid': 33}
	box.writeSAMcatalog(name, mmin=100*box.Melement)
if sf['N_columns'][ii]==41 :
	box.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'c_to_a_500c': 33, 'Ax_500c': 34, 'Ay_500c': 35, 'Az_500c': 36, 'TU': 37, 'M_pe_Behroozi': 38, 'M_pe_Diemer': 39, 'pid': 40}
	box.writeSAMcatalog(name, mmin=100*box.Melement)

print "time used", time.time()-t0, "seconds"

