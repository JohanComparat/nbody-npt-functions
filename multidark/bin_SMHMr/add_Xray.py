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

# open the correct duty cycle tabulate file ...
stellar_mass, duty_cycle_data = n.loadtxt(os.path.join("..", "..", "data", "duty_cycle_0.74230_.txt"), unpack=True)
duty_cycle = interp1d(stellar_mass, duty_cycle_data)

print " set up box, and redshift "
#MD 1 hlist_0.74980_SAM_Nb_0.fits
#MD 25 hlist_0.75440_SAM_Nb_10.fits

def create_catalogs(aexp = 0.74230, env='MD04' , file_type= "hlist"):
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "snapshots", file_type+"_*_stellar_mass.fits" )))
	fileList.sort()
	z = 1./0.74230 -1.
	fileList.sort()
	print fileList
	# set up the stellar mass computation
	sm = StellarMass.StellarMass()
	mhs = n.logspace(10,15,99)

	ratio = sm.SMHMr(mhs,0.)
	stellar_mass = sm.meanSM(mhs,0.)

	# set up the x ray lambda SAR
	logMs = n.arange(6.5,12.5,0.01)
	cdfs_interpolations = []
	#cdfs_interpolations_maxs = []
	XXS = n.arange(32,36.1,0.1)
	for jj,mass in enumerate(logMs):
		norming = xr.Phi_stellar_mass(mass, z)
		cdfs_interpolations.append( interp1d(n.array([xr.Phi_stellar_mass_to_X(X, mass, z) for X in XXS ])/norming, XXS) )
		#cdfs_interpolations_maxs.append( norming )

	cdfs_interpolations = n.array(cdfs_interpolations)
	#cdfs_interpolations_maxs = n.array(cdfs_interpolations_maxs)

	print " loop on the files "
	ii=0
	for fileName in fileList:
		t0=time.time()
		outFile = os.path.join(os.environ[env], "catalogs", os.path.basename(fileName)[:-5] + "_Xray.fits")
		print outFile
		hd = fits.open(fileName)
		
		Nhalo=len(hd[1].data['mvir'])
		randomX = n.random.rand(len(hd[1].data['stellar_mass_Mo13_mvir']))
		indexes = n.searchsorted(logMs,hd[1].data['stellar_mass_Mo13_mvir'])
		lambda_sar_Bo16 = n.array([ cdfs_interpolations[indexes[ii]](randomX[ii]) for ii in range(Nhalo) ])
		
		active_gn = n.array([ duty_cycle(hd[1].data['stellar_mass_Mo13_mvir']) > randomX ])
		
		#Mgal_m200c_Mo13 = norm.rvs( loc = sm.meanSM(10**hd[1].data['mvir'], z), scale = 0.15 )
		#randomY = n.random.rand(len(Mgal_m200c_Mo13))
		#indexesY = n.searchsorted(logMs,Mgal_m200c_Mo13)
		#lambda_sar_Bo16_m200c = n.array([ cdfs_interpolations[indexesY[ii]](randomY[ii]) for ii in range(Nhalo) ])
		#active_gn_m200c = n.array([ cdfs_interpolations_maxs[indexesY[ii]] > randomY[ii] for ii in range(Nhalo) ])
		
		# columns related to Xray AGN
		col01 = fits.Column(name='lambda_sar_Bo16',format='D', array = lambda_sar_Bo16 )
		col02 = fits.Column(name='AGN',format='L', array = active_gn )
		
		#col10 = fits.Column(name='Mgal_m200c_Mo13',format='D', array = Mgal_m200c_Mo13 )
		#col11 = fits.Column(name='lambda_sar_Bo16_M200c',format='D', array = lambda_sar_Bo16_m200c )
		#col12 = fits.Column(name='AGN_M200c',format='L', array = active_gn_m200c )
		
		# columns related to clusters
		col3 = fits.Column(name='Mgas_cluster',format='D', array =n.log10(cl.logM500_to_logMgas(hd[1].data['M500c'], z)))
		col4 = fits.Column(name='kT_cluster',format='D', unit='keV', array =cl.logM500_to_kT(hd[1].data['M500c'], z))
		col5 = fits.Column(name='Lx_bol_cluster',format='D', array =n.log10(cl.logM500_to_L(hd[1].data['M500c'], z)))
		col6 = fits.Column(name='Lx_ce_cluster',format='D', array =n.log10(cl.logM500_to_Lce(hd[1].data['M500c'], z)))

		#define the table hdu 
		colArray = []
		for col in hd[1].columns :
			colArray.append(col)
		
		# AGN Mvir cols
		#colArray.append(col00)
		colArray.append(col01)
		colArray.append(col02)
		# AGN M200c cols
		#colArray.append(col10)
		#colArray.append(col11)
		#colArray.append(col12)
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

create_catalogs(aexp = 0.74230, env='MD04', file_type="hlist")
create_catalogs(aexp = 0.74230, env='MD04', file_type="out")
create_catalogs(aexp = 0.74980, env='MD10', file_type="hlist")
create_catalogs(aexp = 0.75440, env='MD25', file_type="hlist")
create_catalogs(aexp = 0.74980, env='MD10', file_type="out")
create_catalogs(aexp = 0.75440, env='MD25', file_type="out")

