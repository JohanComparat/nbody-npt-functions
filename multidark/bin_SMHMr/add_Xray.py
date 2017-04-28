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

def create_catalogs(env='MD04', volume=400.**3.,  file_type="hlist", aexp='0.74230', out_dir = os.path.join("../../data/")):
	# gets the file list to add the Xray luminosity
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "catalogs", file_type+"*"+aexp+"*stellar_mass.fits" )))
	fileList.sort()
	z = 1./0.74230 -1.
	print fileList
	# opens the duty cycle file_type
	path_to_duty_cycle = os.path.join(out_dir, env+"_"+file_type+"_"+aexp+"_duty_cycle.txt")
	log_stellar_mass, duty_cycle = n.loadtxt(path_to_duty_cycle, unpack="True")
	percentage_active = interp1d(log_stellar_mass, duty_cycle)

	# set up the x ray lambda SAR
	logMs = n.arange(6.5,12.5,0.01)
	cdfs_interpolations = []
	XXS = n.arange(32,36.1,0.1)
	for jj,mass in enumerate(logMs):
		norming = xr.Phi_stellar_mass(mass, z)
		cdfs_interpolations.append( interp1d(n.array([xr.Phi_stellar_mass_to_X(X, mass, z) for X in XXS ])/norming, XXS) )

	cdfs_interpolations = n.array(cdfs_interpolations)

	print " loop on the files "
	for fileName in fileList:
		t0=time.time()
		outFile = os.path.join(os.environ[env], "catalogs", os.path.basename(fileName)[:-5] + "_Xray.fits")
		print outFile
		hd = fits.open(fileName)
		
		stellar_mass = hd[1].data['stellar_mass_Mo13_mvir']
		selection = hd[1].data['stellar_mass_reliable']
		Nhalo=len(stellar_mass)
		
		randomX = n.random.rand(Nhalo)
		active_gn = n.array([ duty_cycle(stellar_mass) > randomX ])
		
		indexes = n.searchsorted(logMs,stellar_mass)
		lambda_sar_Bo16 = n.array([ cdfs_interpolations[indexes[ii]](randomX[ii]) for ii in range(Nhalo) ])
		
		# columns related to Xray AGN
		col1 = fits.Column(name='lambda_sar_Bo16',format='D', array = lambda_sar_Bo16 )
		col2 = fits.Column(name='activity',format='L', array = active_gn )
		
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

create_catalogs(aexp = 0.74230, env='MD04', file_type="out")
create_catalogs(aexp = 0.74980, env='MD10', file_type="out")
create_catalogs(aexp = 0.75440, env='MD25', file_type="out")
#create_catalogs(aexp = 0.74230, env='MD04', file_type="hlist")
#create_catalogs(aexp = 0.74980, env='MD10', file_type="hlist")
#create_catalogs(aexp = 0.75440, env='MD25', file_type="hlist")

