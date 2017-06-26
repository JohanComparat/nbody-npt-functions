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

# dedicated packages
import ClusterScalingRelations
cl = ClusterScalingRelations.ClusterScalingRelations_Mantz2016()
import StellarMass
import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()
dV = -9999

# read the Xray AGN luminosity function and add a condition to reproduce it

def create_catalogs_out(fileList, z, snap_name):
	"""
	Adds Xray emission mass using the Bongiorno et al. 2016 model to the rockstar outputs. 
	"""
	#
	out_duty_cycle = os.path.join(os.environ['MD10'],"duty_cycle", snap_name + "_duty_cycle.txt")
	log_stellar_mass, duty_cycle = n.loadtxt(out_duty_cycle, unpack="True")
	percentage_active = interp1d(n.hstack((-200., 0,n.min(log_stellar_mass)-0.01,log_stellar_mass,n.max(log_stellar_mass)+0.01,15)), n.hstack(( 0., 0., 0., duty_cycle, 0., 0.)))

	# set up the x ray lambda SAR
	logMs = n.arange(4.5,14.5,0.01)
	cdfs_interpolations = []
	XXS = n.arange(32,36.1,0.1)
	for mass in logMs:
		norming = xr.Phi_stellar_mass(mass, z)
		cdfs_interpolations.append( interp1d(n.hstack((n.array([xr.Phi_stellar_mass_to_X(X, mass, z) for X in XXS ])/norming, 1.)), n.hstack((XXS, XXS[-1]+0.1))) )

	cdfs_interpolations = n.array(cdfs_interpolations)
	
	# loops over files
	for fileName in fileList:
		t0=time.time()
		outFile = fileName[:-5]+"_Xray.fits"
		# opens all relevant files
		msFile = fileName[:-5]+"_Ms.fits"
		hd = fits.open(fileName)
		hm = fits.open(msFile)
		
		stellar_mass = hm[1].data['stellar_mass_Mo13_mvir']
		selection = hm[1].data['stellar_mass_reliable']
		Nhalo=len(stellar_mass)
		randomX = n.random.rand(Nhalo)
		active_gn = ( percentage_active(stellar_mass) > randomX )

		indexes = n.searchsorted(logMs,stellar_mass)
		indexes[selection] = n.zeros_like(indexes[selection])
		lambda_sar_Bo16 = n.array([ cdfs_interpolations[indexes[ii]](randomX[ii]) for ii in range(Nhalo) ])

		# columns related to Xray AGN
		col1 = fits.Column(name='lambda_sar_Bo16',format='D', array = lambda_sar_Bo16 )
		col2 = fits.Column(name='activity',format='L', array = active_gn )

		# columns related to clusters
		col3 = fits.Column(name='Mgas_cluster'  ,format='D', array =n.log10(cl.logM500_to_logMgas(hd[1].data['M500c'], z)))
		col4 = fits.Column(name='kT_cluster'    ,format='D', unit='keV', array =cl.logM500_to_kT(hd[1].data['M500c'], z))
		col5 = fits.Column(name='Lx_bol_cluster',format='D', array =n.log10(cl.logM500_to_L(hd[1].data['M500c'], z)))
		col6 = fits.Column(name='Lx_ce_cluster' ,format='D', array =n.log10(cl.logM500_to_Lce(hd[1].data['M500c'], z)))

		#define the table hdu 
		colArray = []
		for col in hd[1].columns :
			colArray.append(col)

		# AGN Mvir cols
		colArray.append(col1)
		colArray.append(col2)
		# Clusters columns
		colArray.append(col3)
		colArray.append(col4)
		colArray.append(col5)
		colArray.append(col6)

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
