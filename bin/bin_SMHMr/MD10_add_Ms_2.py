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

def create_catalogs_out(fileList, z):
	"""
	Adds stellar mass using the Moster et al. 2013 model to the rockstar outputs. 
	"""
	for fileName in fileList:
		t0=time.time()
		outFile = fileName[:-5]+"_Ms.fits"
		hd = fits.open(fileName)
		mean_SM = sm.meanSM(10**hd[1].data['mvir'], z)
		print "mean mgal", mean_SM
		Mgal_mvir_Mo13 = n.array([norm.rvs( loc = mmm, scale = 0.14 ) for mmm in mean_SM ])-n.log10(0.6777)
		print "res  mgal", Mgal_mvir_Mo13
		print "diff mgal - mvir", n.mean(mean_SM-Mgal_mvir_Mo13) 
		print "mean, std magl - mh",n.mean(mean_SM-Mgal_mvir_Mo13), n.std(mean_SM-Mgal_mvir_Mo13)
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
		print time.time()-t0

# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for ii in range(len(summ))[9:18]:
	print summ[ii]
	fileList = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+summ['snap_name'][ii]+'_SAM_Nb_?.fits')))
	#outFile = fileName[:-5]+"_Ms.fits"
	z = summ['redshift'][ii]
	print fileList
	create_catalogs_out(fileList, z)

#create_catalogs(aexp = 0.75440, env='MD25', file_type="hlist")
sys.exit()


if env=='MD04' :
	minMS = 7.2
if env=='MD10' :
	minMS = 9.7
if env=='MD25' :
	minMS = 11.3

