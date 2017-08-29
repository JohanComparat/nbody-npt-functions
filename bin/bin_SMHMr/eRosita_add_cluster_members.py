import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys
import scipy.spatial.ckdtree as t

# specific functions
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d
import random

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
cosmoDS = FlatLambdaCDM(H0=68.46*u.km/u.s/u.Mpc, Om0=0.298734, Ob0=0.046961)
import astropy.constants as constants

path_2_cluster_catalog = os.path.join(os.environ['MD10'], "erosita_4most","MDPL2_FluxProj000_ClustersCombinedModel_withHeader_withExpTime.fits")
hd_cluster = fits.open(path_2_cluster_catalog)

hd_cluster[1].data['snapshot_number']

summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

#second the galaxy member files
for ii, el in enumerate(summ[:-1]):
	snap_N = int(el['snap_name'][:-1])
	bcg = (hd_cluster[1].data['snapshot_number']==snap_N)
	N_clusters = len(bcg.nonzero()[0])
	print snap_N, N_clusters
	#hd_cluster[1].data[cluster_selection]
	file_snaps = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?.fits')))
	file_snaps.sort()
	print file_snaps
	if N_clusters > 0 :
	        file_snaps = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?.fits')))
		file_snaps.sort()
       		print file_snaps
       
		# match all the bcg ids 
		
		# project all the satelittes ... correclty like in the light cone.

		pairs, pairs_N = [], []
		for file_snap in fileList_snaps:
			data = fits.open(file_snap)[1].data
			x, y, z = data[bcg]['x'], data[bcg]['y'], data[bcg]['z']
			radial_distance_2_member = 5. * data[bcg]['rvir'] /1000.
			position_BCG = n.transpose([x, y, z])
			for snap_id, file_snap in enumerate(file_snaps):
				print file_snap
				data = fits.open(file_snap)[1].data
				ok = (data['x']>0)&(data['y']>0)&(data['z']>0)
				treeData = t.cKDTree(n.transpose([ data['x'][ok], data['y'][ok], data['z'][ok] ]),10.0)
				out = n.array([ treeData.query_ball_point(position_BCG[id_bcg], radial_distance_2_member[id_bcg]) for id_bcg in range(len(position_BCG))])
				pairs.append( out )
				pairs_N.append( n.array([len(ppp) for ppp in out]) )
				print len(data['x'][ok]), len(data['x']), time.time(), pairs_N[-1].min(), pairs_N[-1].max()

		# recombine information
		for snap_id, file_snap in enumerate(file_snaps):
			file_out = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(snap_id)+'_4MOST_S5_GAL.fits')
			print( file_out )
			ids_gal = n.hstack((n.array([n.hstack((pairs[snap_id+len(file_snaps)*bcg_counter])).astype('int') for bcg_counter in n.arange(len(fileList_bcg)) ]) ))
			hdu_cols  = fits.ColDefs([
			fits.Column(name='line_number',format='K', array= ids_gal )])
			tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
			prihdr = fits.Header()
			prihdr['HIERARCH nameSnapshot'] = el['snap_name']
			prihdr['batchN'] = snap_id
			prihdr['tracer'] = 'GAL'
			prihdr['author'] = 'JC'
			prihdu = fits.PrimaryHDU(header=prihdr)
			#writes the file
			thdulist = fits.HDUList([prihdu, tb_hdu])
			if os.path.isfile(file_out):
				os.system("rm "+file_out)
			thdulist.writeto(file_out)

