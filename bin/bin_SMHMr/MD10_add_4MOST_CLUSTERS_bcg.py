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

hd = fits.open(os.path.join(os.environ['GIT_NBODY_NPT'], "data/NZ/dndzdom_clusters_Jul17.fits"))
#dn_dv = hd[1].data['DNDZDOM'] / u.deg**2 / cosmoMD.differential_comoving_volume(hd[1].data['Z']).to(u.Mpc**3./u.deg**2)
#*129600. / n.pi
dndv_interpolation = interp1d(n.hstack((-1., 0.,hd[1].data['Z'], 2., 7.)), n.hstack((0., 0.,  hd[1].data['DNDZDOM'], 0., 0.)))
nz_cumulative = n.array([quad(dndv_interpolation, zz, 2.)[0] for zz in zmean]) * 129600. /n.pi
NZ_cumul = interp1d(n.hstack((-1., 0.,zmean, 5., 7.)), n.hstack((0., 0.,  nz_cumulative, 0., 0.)))

#dndv_interpolation_D = interp1d(n.hstack((-1., 0., cosmoMD.comoving_distance(hd[1].data['Z']), cosmoMD.comoving_distance(2.), cosmoMD.comoving_distance(7.))), n.hstack((0., 0., dn_dv, 0., 0.)))

#quad(dndv_interpolation, 0, 2)
#quad(dndv_interpolation_D, 0, cosmoMD.comoving_distance(2.).value)

#zmin, zmax, zmean, dN_dV_bcg_data, dN_dV_gal_data = n.loadtxt(os.path.join(os.environ['GIT_NBODY_NPT'], "data/NZ/nz_5.txt"), unpack=True)

#NZ_cumul = lambda z : 10**(5. - 1.9*z)

dVolume = (cosmoMD.comoving_volume(zmax) - cosmoMD.comoving_volume(zmin))
dN_dV_bcg_data = (NZ_cumul(zmin)-NZ_cumul(zmax))/dVolume
dN_dV_bcg = interp1d(n.hstack((-1., zmin[0], zmean, zmax[-1], 3.)), n.hstack((0., 0., dN_dV_bcg_data, 0., 0.)) )


N_bcg_per_snap = (volume * dN_dV_bcg(z_snap) ).astype('int')
#N_bcg_per_snap2 = (volume * dndv_interpolation(z_snap) ).astype('int')
# first the BCG files
L_min = n.zeros(len(summ))
for ii in n.arange(len(summ)):
	N_cluster = N_bcg_per_snap[ii]
	el = summ[ii]
	print(ii, el, N_cluster)
	fileList_snap_X = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?_Xray.fits')))
	fileList_snap_X.sort()
	if N_cluster > 1 :
		# LRG: select on Ms
		LX = n.array([fits.open(file)[1].data['Lx_bol_cluster'] for file in fileList_snap_X])

		all_LX = n.hstack((LX))
		all_LX_sort_id = n.argsort(all_LX)
		min_LX = n.min(all_LX[all_LX_sort_id[-N_cluster-1:]])
		L_min[ii] = min_LX
		print(min_LX)
		id_superset = n.where(LX>min_LX)

		for num in list(set(id_superset[0])):
			file_out = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(num)+'_4MOST_S5_BCG.fits')
			bcgs = id_superset[1][ (id_superset[0]==num) ]
			print( file_out, bcgs)
			hdu_cols  = fits.ColDefs([
			fits.Column(name='line_number',format='K', array= bcgs )])
			tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
			prihdr = fits.Header()
			prihdr['HIERARCH nameSnapshot'] = el['snap_name']
			prihdr['batchN'] = num
			prihdr['tracer'] = 'BCG'
			prihdr['author'] = 'JC'
			prihdu = fits.PrimaryHDU(header=prihdr)
			#writes the file
			thdulist = fits.HDUList([prihdu, tb_hdu])
			if os.path.isfile(file_out):
				os.system("rm "+file_out)
			thdulist.writeto(file_out)
