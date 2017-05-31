# overall python packages
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys
t0=time.time()
		

def create_AGN_catalog(env='MD10', file_type="out", aexp='0.74230'):
	# gets the file list to add the Xray luminosity
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", file_type+"_"+aexp+"_*stellar_mass_Xray.fits" )))
	fileList.sort()
	z = 1./0.74230 -1.
	print fileList
	print " loop on the files "
	outFile = os.path.join(os.environ[env], "catalogs", file_type + "_"+ aexp + "_AGN.fits")
	print "out in ", outFile
	print fileList[0]
	hd = fits.open(fileList[0])
	agn = (hd[1].data['activity'] == True)
	data_all = hd[1].data[agn]
	for fileName in fileList[1:]:
		print fileName
		hd = fits.open(fileName)
		agn = (hd[1].data['activity'] == True)
		data_all = n.hstack((data_all, hd[1].data[agn]))
		
		tb_hdu = fits.BinTableHDU.from_columns( data_all )
		#define the header
		prihdr = fits.Header()
		prihdr['author'] = 'JC'
		prihdr['redshift'] = z
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		if os.path.isfile(outFile):
			os.system("rm "+outFile)
			
		thdulist.writeto(outFile)
		

create_AGN_catalog(env='MD10', file_type="out"  , aexp='0.74980')
print time.time()-t0, "seconds"
