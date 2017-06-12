# overall python packages
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n

# specific functions
from scipy.stats import norm

# dedicated packages
import StellarMass
sm = StellarMass.StellarMass()

print " set up box, and redshift "
#MD 1 hlist_0.74980_SAM_Nb_0.fits
#MD 25 hlist_0.75440_SAM_Nb_10.fits


def create_catalogs(aexp = 0.74230, env='MD04' , file_type= "hlist", dV=-99.99):
	if env=='MD04' :
		minMS = 7.2
	if env=='MD10' :
		minMS = 9.7
	if env=='MD25' :
		minMS = 11.3

	fileList = n.array(glob.glob(os.path.join(os.environ[env], "snapshots", file_type+"_*_SAM_Nb_*.fits" )))
	fileList.sort()
	z = 1./0.74230 -1.
	fileList.sort()
	print fileList
	for fileName in fileList:
		t0=time.time()
		outFile = os.path.join(os.environ[env], "catalogs", os.path.basename(fileName)[:-5] + "_stellar_mass.fits")
		print outFile
		hd = fits.open(fileName)
		
		Mgal_mvir_Mo13 = norm.rvs( loc = sm.meanSM(10**hd[1].data['mvir'], z), scale = 0.15 )-n.log10(0.6777)
		sel = (Mgal_mvir_Mo13>minMS)
		
		col00 = fits.Column(name='stellar_mass_Mo13_mvir',format='D', unit='logMsun', array = Mgal_mvir_Mo13 )
		col01 = fits.Column(name='stellar_mass_reliable', format='L', array = sel )
		
		#define the table hdu 
		colArray = []
		for col in hd[1].columns :
			colArray.append(col)
		
		# Mvir stellar mass
		colArray.append(col00)
		colArray.append(col01)
		
		hdu_cols  = fits.ColDefs(colArray)
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		
		#define the header
		prihdr = fits.Header()
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		if os.path.isfile(outFile):
			os.system("rm "+outFile)
			
		thdulist.writeto(outFile)
		print time.time()-t0

create_catalogs(aexp = 0.74230, env='MD04', file_type="hlist")
create_catalogs(aexp = 0.74230, env='MD04', file_type="out")
create_catalogs(aexp = 0.74980, env='MD10', file_type="hlist")
create_catalogs(aexp = 0.75440, env='MD25', file_type="hlist")
create_catalogs(aexp = 0.74980, env='MD10', file_type="out")
create_catalogs(aexp = 0.75440, env='MD25', file_type="out")

