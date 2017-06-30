import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

# specific functions
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d
import random

summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	
volume = (1000/0.6777)**3.

zmin, zmax, zmean, dN_dV_lrg_data, dN_dV_elg_data, dN_dV_qso_data = n.loadtxt(os.path.join(os.environ['GIT_NBODY_NPT'], "data/NZ/nz_8.txt"), unpack=True)

dN_dV_elg = interp1d(n.hstack((-1., zmin[0], zmean, zmax[-1], 3.)), n.hstack((0., 0., dN_dV_elg_data, 0., 0.)) )
dN_dV_lrg = interp1d(n.hstack((-1., zmin[0], zmean, zmax[-1], 3.)), n.hstack((0., 0., dN_dV_lrg_data, 0., 0.)) )
dN_dV_qso = interp1d(n.hstack((-1., zmin[0], zmean, zmax[-1], 3.)), n.hstack((0., 0., dN_dV_qso_data, 0., 0.)) )

# number of tracer per snapshot (complete volume)
z_snap = summ['redshift']

N_elg_per_snap = (volume * dN_dV_elg(z_snap) ).astype('int')
N_lrg_per_snap = (volume * dN_dV_lrg(z_snap) ).astype('int')
N_qso_per_snap = (volume * dN_dV_qso(z_snap) ).astype('int')

# determine dz
#z_middle = (summ['redshift'][1:]+summ['redshift'][:-1])*0.5
#z_mins = n.hstack((summ['redshift'][0], z_middle))
#z_maxs = n.hstack((z_middle, summ['redshift'][-1]))
#dz = z_maxs - z_mins

for ii, el in enumerate(summ):
	print(N_lrg_per_snap[ii], N_elg_per_snap[ii], N_qso_per_snap[ii])
	N_lrg, N_elg, N_qso = N_lrg_per_snap[ii], N_elg_per_snap[ii], N_qso_per_snap[ii]

	fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?.fits')))
	fileList_snap.sort()
	fileList_snap_MS = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?_Ms.fits')))
	fileList_snap_MS.sort()

	if N_lrg > 10 :
		# LRG: select on Ms
		MS = n.array([fits.open(file)[1].data['stellar_mass_Mo13_mvir'] for file in fileList_snap_MS])

		all_MS = n.hstack((MS))
		all_MS_sort_id = n.argsort(all_MS)
		min_mass = all_MS[all_MS_sort_id[-N_lrg*2-1]]

		id_superset = n.where(MS>min_mass)
		rds = n.random.rand(len(id_superset[0]))

		for num in list(set(id_superset[0])):
			file_out = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(num)+'_4MOST_S8_LRG.fits')
			lrgs = id_superset[1][ (id_superset[0]==num) & (rds>0.5) ]
			print( file_out, lrgs)
			hdu_cols  = fits.ColDefs([
			fits.Column(name='line_number',format='I', array= lrgs )])
			tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
			prihdr = fits.Header()
			prihdr['HIERARCH nameSnapshot'] = el['snap_name']
			prihdr['batchN'] = num
			prihdr['tracer'] = 'LRG'
			prihdr['author'] = 'JC'
			prihdu = fits.PrimaryHDU(header=prihdr)
			#writes the file
			thdulist = fits.HDUList([prihdu, tb_hdu])
			if os.path.isfile(file_out):
				os.system("rm "+file_out)
			thdulist.writeto(file_out)

	if N_elg > 10 or N_qso > 10  :
		MH = n.array([fits.open(file)[1].data['mvir'] for file in fileList_snap])

		if N_elg > 10 :
			# ELG select on Mvir
			#p_elg=[10**(12.2),0.25]
			mh_mean, mh_scatter = 12.2, 0.25
			mh_bins = n.arange(mh_mean -2*mh_scatter, mh_mean +2*mh_scatter+0.1, 0.1)
			mh_bins_pos = 0.5*(mh_bins[1:]+mh_bins[:-1])
			proba = lambda x : norm.pdf(x, loc=12.2,scale=0.25)
			proba_norm = proba(mh_bins_pos).sum()
			N_2_select_per_bin = (N_elg*proba(mh_bins_pos)/proba_norm).astype('int')

			id_0 = []
			id_1 = []
			for id_bin in range(len(mh_bins)-1):
				id_superset = n.where( (MH > mh_bins[id_bin]) &( MH < mh_bins[id_bin+1]) )
				N_avail = id_superset[0].shape[0]
				rds = n.random.rand(len(id_superset[0]))
				bin_selection = (rds < N_2_select_per_bin[id_bin]*1./N_avail)
				id_0.append(id_superset[0][bin_selection])
				id_1.append(id_superset[1][bin_selection])

			id_0 = n.hstack((id_0))
			id_1 = n.hstack((id_1))

			for num in list(set(id_0)):
				file_out = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(num)+'_4MOST_S8_ELG.fits')
				elgs = id_1[ (id_0==num) ]
				print( file_out, elgs)
				hdu_cols  = fits.ColDefs([
				fits.Column(name='line_number',format='I', array= elgs )])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = el['snap_name']
				prihdr['batchN'] = num
				prihdr['tracer'] = 'ELG'
				prihdr['author'] = 'JC'
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				if os.path.isfile(file_out):
					os.system("rm "+file_out)
				thdulist.writeto(file_out)

		if N_qso > 10 :
			# QSO select on Mvir
			#p_qso=[10**(12.7),0.25]
			mh_mean, mh_scatter = 12.7, 0.25
			mh_bins = n.arange(mh_mean -2*mh_scatter, mh_mean +2*mh_scatter+0.1, 0.1)
			mh_bins_pos = 0.5*(mh_bins[1:]+mh_bins[:-1])
			proba = lambda x : norm.pdf(x, loc=12.2,scale=0.25)
			proba_norm = proba(mh_bins_pos).sum()
			N_2_select_per_bin = (N_qso*proba(mh_bins_pos)/proba_norm).astype('int')

			id_0 = []
			id_1 = []
			for id_bin in range(len(mh_bins)-1):
				id_superset = n.where( (MH > mh_bins[id_bin]) &( MH < mh_bins[id_bin+1]) )
				N_avail = id_superset[0].shape[0]
				rds = n.random.rand(len(id_superset[0]))
				bin_selection = (rds < N_2_select_per_bin[id_bin]*1./N_avail)
				id_0.append(id_superset[0][bin_selection])
				id_1.append(id_superset[1][bin_selection])

			id_0 = n.hstack((id_0))
			id_1 = n.hstack((id_1))

			for num in list(set(id_0)):
				file_out = os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_'+str(num)+'_4MOST_S8_QSO.fits')
				qsos = id_1[ (id_0==num) ]
				print( file_out, qsos)
				hdu_cols  = fits.ColDefs([
				fits.Column(name='line_number',format='I', array= qsos )])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = el['snap_name']
				prihdr['batchN'] = num
				prihdr['tracer'] = 'QSO'
				prihdr['author'] = 'JC'
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				if os.path.isfile(file_out):
					os.system("rm "+file_out)
				thdulist.writeto(file_out)


