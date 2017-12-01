# overall python packages
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys
from scipy.interpolate import interp1d
# read the Xray AGN luminosity function and add a condition to reproduce it
# obscuration law
obscuration_z_grid, obscuration_nh_grid, obscuration_fraction_obs_erosita = n.loadtxt( os.path.join( os.environ['GIT_NBODY_NPT'], "data", "AGN", "fraction_observed_by_erosita_due_2_obscuration.txt"), unpack=True)
nh_vals = 10**n.arange(-2,4,0.05)
z_vals = 10**n.arange(-3,0.68,0.025)
obscuration_interpolation_grid = n.array([ 
  interp1d(
    n.hstack((obscuration_nh_grid[ (obscuration_z_grid==zz) ], 26.)), 
    n.hstack((obscuration_fraction_obs_erosita[( obscuration_z_grid==zz) ], obscuration_fraction_obs_erosita[( obscuration_z_grid==zz) ][-1]))
	      ) 
  for zz in z_vals])

def create_catalog(snap_name, z):
	"""
	Creates summary catalog for AGN only"""
	#
	fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+snap_name+'_SAM_Nb_?.fits')))
	fileList_ms = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+snap_name+'_SAM_Nb_?_Ms.fits')))
	fileList_Xray = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+snap_name+'_SAM_Nb_?_LSAR.fits')))
	fileList_DC = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+snap_name+'_SAM_Nb_?_DC.fits')))
	
	fileList_snap.sort()
	fileList_ms.sort()
	fileList_Xray.sort()
	fileList_DC.sort()

	out_snap = os.path.join(os.environ['MD10'], "catalogs", 'out_'+snap_name + "_AGN_snapshot.fits")
	out_ms = os.path.join(os.environ['MD10'], "catalogs", 'out_'+snap_name + "_AGN_Ms.fits")
	out_xray = os.path.join(os.environ['MD10'], "catalogs", 'out_'+snap_name + "_AGN_LSAR.fits")
	#out_dc = os.path.join(os.environ['MD10'], "catalogs", 'out_'+snap_name + "_AGN_DC.fits")

	# loops over files
	dat_snap = []
	dat_ms = []
	dat_xray = []
	LX = []
	index = n.searchsorted(z_vals, z)

	for fileSnap, fileMs, fileXray, fileDC in zip(fileList_snap, fileList_ms, fileList_Xray, fileList_DC):
		print fileSnap
		print fileMs
		print fileXray
		print fileDC
		
		hd = fits.open(fileSnap)
		hm = fits.open(fileMs)
		hx = fits.open(fileXray)
		hc = fits.open(fileDC)
		
		MS = hm[1].data['stellar_mass_Mo13_mvir']
		SAR = hx[1].data['lambda_sar_Bo16']
		agn = (hc[1].data['activity']==1)&(MS>0)
		#lognh = hx[1].data['log_NH_Buchner2017']

		dat_snap.append( hd[1].data[agn])
		dat_ms.append( hm[1].data[agn])
		dat_xray.append( hx[1].data[agn])
		
		#percent_observed = obscuration_interpolation_grid[index]( lognh ) 
		LX.append( MS[agn] + SAR[agn] )



	hdu_cols  = fits.ColDefs(n.hstack((dat_snap)))
	print "snap", n.hstack((dat_snap)).shape
	tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
	#define the header
	prihdr = fits.Header()
	prihdr['author'] = 'JC'
	prihdr['info'] = 'snapshot'
	prihdu = fits.PrimaryHDU(header=prihdr)
	#writes the file
	thdulist = fits.HDUList([prihdu, tb_hdu])
	if os.path.isfile(out_snap):
		os.system("rm "+out_snap)
		
	thdulist.writeto(out_snap)

	hdu_cols  = fits.ColDefs(n.hstack((dat_ms)))
	print "ms", n.hstack((dat_ms)).shape
	tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
	#define the header
	prihdr = fits.Header()
	prihdr['author'] = 'JC'
	prihdr['info'] = 'ms'
	prihdu = fits.PrimaryHDU(header=prihdr)
	#writes the file
	thdulist = fits.HDUList([prihdu, tb_hdu])
	if os.path.isfile(out_ms):
		os.system("rm "+out_ms)
		
	thdulist.writeto(out_ms)

	hdu_cols  = fits.ColDefs(n.hstack((dat_xray)))
	hdu_cols.add_col( fits.Column(name='LX_2_10_keV',format='D', array= n.hstack((LX)) ))
	print "xray", n.hstack((dat_xray)).shape
	tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
	#define the header
	prihdr = fits.Header()
	prihdr['author'] = 'JC'
	prihdr['info'] = 'xray'
	prihdu = fits.PrimaryHDU(header=prihdr)
	#writes the file
	thdulist = fits.HDUList([prihdu, tb_hdu])
	if os.path.isfile(out_xray):
		os.system("rm "+out_xray)
		
	thdulist.writeto(out_xray)



# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for el in summ[:1]:
	print el
	create_catalog(el['snap_name'], el['redshift'])
