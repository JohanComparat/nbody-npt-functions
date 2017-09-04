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

from numpy import exp, random

logmbh = lambda logm200c : 8.18 + 1.55 * ( logm200c - 13. )

LX_Z_lo, FRACTION_Z_lo = n.loadtxt(os.path.join(os.environ['GIT_NBODY_NPT'],'data','fraction_zgt0_zlt10.txt' ), unpack=True)
LX_Z_me, FRACTION_Z_me = n.loadtxt(os.path.join(os.environ['GIT_NBODY_NPT'],'data','fraction_zgt10_zlt25.txt'), unpack=True)
LX_Z_hi, FRACTION_Z_hi = n.loadtxt(os.path.join(os.environ['GIT_NBODY_NPT'],'data','fraction_zgt25_zlt50.txt'), unpack=True)
zmin_obscur = n.array([-0.01,1.0,2.5])
zmax_obscur = n.array([1.0,2.5,10.0])

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
	
	if z > zmin_obscur[0] and z<=zmax_obscur[0]:
		obscured_fraction_interpolated = interp1d(n.hstack(( 0, LX_Z_lo, 60 )), n.hstack(( FRACTION_Z_lo[0],FRACTION_Z_lo,FRACTION_Z_lo[-1] )) )
	elif z > zmin_obscur[1] and z<=zmax_obscur[1]:
		obscured_fraction_interpolated = interp1d(n.hstack(( 0, LX_Z_me, 60 )), n.hstack(( FRACTION_Z_me[0],FRACTION_Z_me,FRACTION_Z_me[-1] )) )
	else:
		obscured_fraction_interpolated = interp1d(n.hstack(( 0, LX_Z_hi, 60 )),  n.hstack(( FRACTION_Z_hi[0],FRACTION_Z_hi,FRACTION_Z_hi[-1] )) )
		

	# loops over files
	for fileName in fileList:
		t0=time.time()
		outFile = fileName[:-5]+"_Xray.fits"
		# opens all relevant files
		msFile = fileName[:-5]+"_Ms.fits"
		hd = fits.open(fileName)
		hm = fits.open(msFile)
		
		stellar_mass = hm[1].data['stellar_mass_Mo13_mvir']
		selection = (hm[1].data['stellar_mass_Mo13_mvir']>0) # hm[1].data['stellar_mass_reliable']
		Nhalo=len(stellar_mass)
		randomX = n.random.rand(Nhalo)
		active_gn = ( percentage_active(stellar_mass) > randomX )
		# lambda SAR addition
		indexes = n.searchsorted(logMs,stellar_mass)
		indexes[selection] = n.zeros_like(indexes[selection])
		lambda_sar_Bo16 = n.array([ cdfs_interpolations[indexes[ii]](randomX[ii]) for ii in range(Nhalo) ])
		# obscuration Merloni et al. 2015 as a function of luminosity
		# randomObscuration = n.random.rand(Nhalo)		
		# obscured = (xr.obscured_fraction_optical_Merloni2015(lambda_sar_Bo16 + stellar_mass) < randomObscuration )
		
		# obscuration Buchner et al. 2015 + 2017
		# add the log NH of the logNH_host
		# 35 % have a thick obscuration 24 - 26
		# 65 % have a thin obscuration that depends on stellar mass and Xray luminosity
		logNH = random.uniform(20, 22,Nhalo)
		obs_type = n.zeros(Nhalo)
		# 35% of thick, 24-26
		randomNH = n.random.rand(Nhalo)
		thick_obscuration = (randomNH < 0.35)
		thin_obscuration = (randomNH >= 0.35)
		logNH[thick_obscuration] = random.uniform(24, 26, len(logNH[thick_obscuration]))
		obs_type[thick_obscuration] = n.ones_like(logNH[thick_obscuration])*2
		# the thin : about 40 % are thin whatever happens: 22-24
		logNH_host_mean = 21.7 + (stellar_mass - 9.5)*0.38
		logNH_host = random.normal(logNH_host_mean, 0.5)
		logNH[(thin_obscuration)&(logNH_host>=22)] = random.uniform(22, 24, len(logNH[(thin_obscuration)&(logNH_host>=22)]))
		obs_type[(thin_obscuration)&(logNH_host>=22)] =  n.ones_like(logNH[(thin_obscuration)&(logNH_host>=22)])
		# a few more are thin depending on their Xray luminosity: 22-24
		rest = (thin_obscuration)&(logNH_host<22)
		randomNH2 = n.random.rand(Nhalo)
		rest_obscured = (rest)&(randomNH2 < obscured_fraction_interpolated(lambda_sar_Bo16 + stellar_mass))
		logNH[(rest_obscured)] = random.uniform(22, 24, len(logNH[(rest_obscured)]))
		obs_type[(rest_obscured)] =  n.ones_like(logNH[(rest_obscured)])
		# the rest has: 20-22 by default
		print logNH, obs_type
		
		# columns related to Xray AGN
		col1 = fits.Column(name='lambda_sar_Bo16',format='D', array = lambda_sar_Bo16 )
		col2 = fits.Column(name='activity',format='L', array = active_gn )
		col3 = fits.Column(name='obscuration_type_Buchner2017',format='K', array = obs_type.astype('int') )
		#col4 = fits.Column(name='M_BH'  ,format='D', array =logmbh(hd[1].data['M200c'])))
		col4 = fits.Column(name='log_NH_Buchner2017'  ,format='D', array = logNH)
		

		# columns related to clusters
		col5 = fits.Column(name='Mgas_cluster'  ,format='D', array =n.log10(cl.logM500_to_logMgas(hd[1].data['M500c'], z)))
		col6 = fits.Column(name='kT_cluster'    ,format='D', unit='keV', array =cl.logM500_to_kT(hd[1].data['M500c'], z))
		col7 = fits.Column(name='Lx_bol_cluster',format='D', array =n.log10(cl.logM500_to_L(hd[1].data['M500c'], z)))
		col8 = fits.Column(name='Lx_ce_cluster' ,format='D', array =n.log10(cl.logM500_to_Lce(hd[1].data['M500c'], z)))

		#define the table hdu 
		colArray = [col1]
		#for col in hd[1].columns :
			#colArray.append(col)

		# AGN Mvir cols
		colArray.append(col2)
		colArray.append(col3)
		# Clusters columns
		colArray.append(col4)
		colArray.append(col5)
		colArray.append(col6)
		colArray.append(col7)
		colArray.append(col8)

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

for el in summ[27:36]:
	print el
	fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?.fits')))
	fileList_snap.sort()
	print fileList_snap
	create_catalogs_out(fileList_snap, el['redshift'], el['snap_name'])
