# overall python packages
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys 
# specific functions
from scipy.stats import norm
# dedicated packages
import StellarMass
sm = StellarMass.StellarMass()

def create_catalogs_out(fileList, z, minMS = 8.7):
	"""
	Adds stellar mass using the Moster et al. 2013 model to the rockstar outputs. 
	"""
	for fileName in fileList:
		t0=time.time()
		outFile = fileName[:-5]+"_Ms.fits"
		hd = fits.open(fileName)
		Mgal_mvir_Mo13 = norm.rvs( loc = sm.meanSM(10**hd[1].data['mvir'], z), scale = 0.15 )-n.log10(0.6777)
		sel = (Mgal_mvir_Mo13>minMS)
		
		col00 = fits.Column(name='stellar_mass_Mo13_mvir',format='D', unit='logMsun', array = Mgal_mvir_Mo13 )
		col01 = fits.Column(name='stellar_mass_reliable', format='L', array = sel )

		#define the table hdu 
		colArray = []
		colArray.append(hd[1].columns[0])
		# Mvir stellar mass
		colArray.append(col00)
		colArray.append(col01)

		hdu_cols  = fits.ColDefs(colArray)
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )

		#define the header
		prihdr = fits.Header()
		prihdr['author'] = 'JC'
		prihdr['SAMfile'] = os.path.basename(fileName)
		prihdr['minMS'] = minMS
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		if os.path.isfile(outFile):
			os.system("rm "+outFile)

		thdulist.writeto(outFile)
		print time.time()-t0

# open the output file_type
summary_file_name = os.path.join(os.environ["NUGC_DIR"], 'output_NUGC.fits')
summ = fits.open(summary_file_name)[1].data	

for ii in range(len(summ)):
	print summ[ii]
	fileList = n.hstack((n.array(glob.glob(os.path.join(os.environ["NUGC_DIR"], 'work_agn', summ['snap_name'][ii][:-5]+'_SAM_Nb_?.fits'))), n.array(glob.glob(os.path.join(os.environ["NUGC_DIR"], 'work_agn', summ['snap_name'][ii][:-5]+'_SAM_Nb_??.fits'))) ))
	z = summ['redshift'][ii]
	print fileList
	create_catalogs_out(fileList, z)

