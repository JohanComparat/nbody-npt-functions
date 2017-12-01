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
#import StellarMass

meanSM= lambda Mh, z : n.log10(Mh * 2. * ( 0.0351 - 0.0247 * z/(1.+z)) / ((Mh/ (10**(11.79 + 1.5 * z/(1.+z))) )**(- 0.9 + 0.5  * z/(1.+z)) + ( Mh /(10**(11.79 + 1.5 * z/(1.+z))) )**(0.67 + 0.2 * z/(1.+z)) ) )

fun = lambda mmm : norm.rvs( loc = mmm, scale = 0.15 )


def create_catalogs_out(fileList, z):
	"""
	Adds stellar mass using the Moster et al. 2013 model to the rockstar outputs. 
	"""
	for fileName in fileList:
		t0=time.time()
		outFile = fileName[:-5]+"_Ms.fits"
		hd = fits.open(fileName)
		mean_SM = meanSM(10**hd[1].data['mvir']/0.6777, z)
		#print "mean mgal", mean_SM
		Mgal_mvir_Mo13 = n.array([fun(el) for el in mean_SM]) # n.array(pool.starmap( fun, mean_SM ))
		#print "res  mgal", Mgal_mvir_Mo13
		#print "diff mgal - mvir", n.mean(mean_SM-Mgal_mvir_Mo13) 
		#print "mean, std magl - mh",n.mean(mean_SM-Mgal_mvir_Mo13), n.std(mean_SM-Mgal_mvir_Mo13)
		sel = (hd[1].data['mvir']>0)
		
		Mgal_mvir_Mo13[sel==False] = n.zeros_like(Mgal_mvir_Mo13[sel==False])
		
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
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		if os.path.isfile(outFile):
			os.system("rm "+outFile)

		thdulist.writeto(outFile)
		print( time.time()-t0)

# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for ii in range(len(summ))[18:27]:
	print( summ[ii])
	fileList = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+summ['snap_name'][ii]+'_SAM_Nb_?.fits')))
	#outFile = fileName[:-5]+"_Ms.fits"
	z = summ['redshift'][ii]
	print( fileList)
	create_catalogs_out(fileList, z)


