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
	out_duty_cycle = os.path.join(os.environ['MD10'],"duty_cycle", "out_" + snap_name + "_duty_cycle.txt")
	log_stellar_mass, duty_cycle = n.loadtxt(out_duty_cycle, unpack="True")
	percentage_active = interp1d(n.hstack((-200., 0,n.min(log_stellar_mass)-0.01,log_stellar_mass,n.max(log_stellar_mass)+0.01,15)), n.hstack(( 0., 0., 0., duty_cycle, 0., 0.)))
	
	# loops over files
	for fileName in fileList:
		t0=time.time()
		outFile = fileName[:-5]+"_DC.fits"
		# opens all relevant files
		msFile = fileName[:-5]+"_Ms.fits"
		hd = fits.open(fileName)
		hm = fits.open(msFile)
		
		logM = hm[1].data['stellar_mass_Mo13_mvir']
		agn_random_number = n.random.random(len(logM))
		activity = n.zeros(len(logM))
		proba = percentage_active(logM)
		activity[agn_random_number < proba] = n.ones_like(activity[agn_random_number < proba])

		
		# columns related to Xray AGN
		col1 = fits.Column(name='activity',format='K', array = activity )
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

for el in summ:
	print el
	fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?.fits')))
	fileList_snap.sort()
	print fileList_snap
	create_catalogs_out(fileList_snap, el['redshift'], el['snap_name'])
