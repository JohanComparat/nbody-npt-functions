# overall python packages
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

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)



# read the Xray AGN luminosity function and add a condition to reproduce it

def create_catalogs_out(fileList, z, snap_name):
	"""
	Adds Xray emission mass using the Bongiorno et al. 2016 model to the rockstar outputs. 
	"""
	def f_lambda_sar( DATA ):
	  logM, log_lambda_SAR = DATA
	  log_lambda_SAR_var = 10**( log_lambda_SAR - 33.8 + 0.48 * (logM - 11.) )
	  return 1. / ( log_lambda_SAR_var**(1.01 - 0.58 * (z - 1.1)) + log_lambda_SAR_var**(3.72) )

	dl=0.01
	log_lambda_SAR_values = n.arange(32-dl,36+2*dl,dl)
	
	# loops over files
	for fileName in fileList:
		t0=time.time()
		outFile = fileName[:-5]+"_LSAR.fits"
		# opens all relevant files
		msFile = fileName[:-5]+"_Ms.fits"
		hd = fits.open(fileName)
		hm = fits.open(msFile)
		
		logM = hm[1].data['stellar_mass_Mo13_mvir']
		agn_random_number = n.random.random(len(logM))
		log_lSAR = n.zeros(len(logM))

		t0 = time.time()
		ii0=0
		ii_step=12000
		for ii0 in n.arange(0, len(logM), ii_step):
			ii1=ii0+ii_step
			X,Y = n.meshgrid(logM[ii0:ii1], log_lambda_SAR_values)
			#Z = n.ones_like(X)*z
			probas_un = f_lambda_sar([ X, Y])#, Z ])
			norm = n.sum(probas_un, axis=0)
			probas = probas_un / norm
			cmat = n.array([ agn_random_number[ii0:ii1] > n.sum(probas.T[:,jj:], axis=1) for jj in n.arange(len(log_lambda_SAR_values)) ])
			#print(cmat.shape)#, cmat[0])
			#print(cmat.T[1])
			print(ii0, len(logM), cmat.shape, time.time()-t0)
			values = log_lambda_SAR_values[n.array([n.min(n.where(cmat.T[jj]==True)) for jj in n.arange(len(cmat.T)) ])]
			#print(values.shape, values[:10])
			log_lSAR[ii0:ii1] = values
			
		
		# columns related to Xray AGN
		col1 = fits.Column(name='lambda_sar_Bo16',format='D', array = log_lSAR )
		col1b = fits.Column(name='agn_random_number',format='D', array = agn_random_number )
		
		#define the table hdu 
		colArray = [col1]
		colArray.append(col1b)
		#for col in hd[1].columns :
			#colArray.append(col)
			
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


# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for el in summ[27:36]:#[27:36]:
	print el
	fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?.fits')))
	fileList_snap.sort()
	print fileList_snap
	create_catalogs_out(fileList_snap, el['redshift'], el['snap_name'])
