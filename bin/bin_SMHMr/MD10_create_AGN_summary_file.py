# overall python packages
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

# read the Xray AGN luminosity function and add a condition to reproduce it

def create_catalog(snap_name, z):
	"""
	Creates summary catalog for AGN only"""
	#
	fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+snap_name+'_SAM_Nb_?.fits')))
	fileList_ms = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+snap_name+'_SAM_Nb_?_Ms.fits')))
	fileList_Xray = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+snap_name+'_SAM_Nb_?_Xray.fits')))
	fileList_snap.sort()
	fileList_ms.sort()
	fileList_Xray.sort()
	
	out_snap = os.path.join(os.environ['MD10'], "catalogs", 'out_'+snap_name + "_AGN_snapshot.fits")
	out_ms = os.path.join(os.environ['MD10'], "catalogs", 'out_'+snap_name + "_AGN_Ms.fits")
	out_xray = os.path.join(os.environ['MD10'], "catalogs", 'out_'+snap_name + "_AGN_Xray.fits")
	
	# loops over files
	dat_snap = []
	dat_ms = []
	dat_xray = []
	
	for fileSnap, fileMs, fileXray in zip(fileList_snap, fileList_ms, fileList_Xray):
		print fileSnap
		print fileMs
		print fileXray
		hd = fits.open(fileSnap)
		hm = fits.open(fileMs)
		hx = fits.open(fileXray)
		
		agn = (hx[1].data['activity'])
		
		dat_snap.append( hd[1].data[agn])
		dat_ms.append( hm[1].data[agn])
		dat_xray.append( hx[1].data[agn])


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

for el in summ:
	print el
	create_catalog(el['snap_name'], el['redshift'])
