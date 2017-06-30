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

dz = 0.05
bins = n.arange(0,5,dz)
zmin, zmax, zmean = bins[:-1], bins[1:], (bins[:-1]+bins[1:])/2.

summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	
volume = (1000/0.6777)**3.
z_snap = summ['redshift']

zmin, zmax, zmean, dN_dV_bcg_data, dN_dV_gal_data = n.loadtxt(os.path.join(os.environ['GIT_NBODY_NPT'], "data/NZ/nz_5.txt"), unpack=True)

NZ_cumul = lambda z : 10**(5.-2.*z)

dVolume = (cosmoMD.comoving_volume(zmax) - cosmoMD.comoving_volume(zmin))
dN_dV_bcg_data = (NZ_cumul(zmin)-NZ_cumul(zmax))/dVolume
dN_dV_bcg = interp1d(n.hstack((-1., zmin[0], zmean, zmax[-1], 3.)), n.hstack((0., 0., dN_dV_bcg_data, 0., 0.)) )

N_bcg_per_snap = (volume * dN_dV_bcg(z_snap) ).astype('int')

#second the galaxy member files
for ii, el in enumerate(summ):
	fileList_bcg = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?_4MOST_S5_BCG.fits')))
	fileList_bcg.sort()
	file_snaps = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?.fits')))
	file_snaps.sort()
	print fileList_bcg
	if len(fileList_bcg)>0 :
		pairs, pairs_N = [], []

		for file_bcg in fileList_bcg:
			print file_bcg
			bcg = fits.open(file_bcg)[1].data['line_number']
			num = os.path.basename(file_bcg).split('_')[-4]
			file_snap = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(num)+'.fits')
			data = fits.open(file_snap)[1].data
			x, y, z = data[bcg]['x'], data[bcg]['y'], data[bcg]['z']
			radial_distance_2_member = 3. * data[bcg]['rvir'] /1000.
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

