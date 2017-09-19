
"""
.. class:: ClusterScalingRelations

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class clusterScalingRelations is a wrapper to add cluster physics to the Multidark simulations results / outputs.

References
----------

 * M10: Mantz et al. 2010b
 * Z14: Zandanel et al. 2014
 * M16: Mantz et al. arxiv. 1606.03407 Weighing the Giants V: Galaxy Cluster Scaling Relations
 * M15: Mantz et al. arxiv. 1509.01322 Cosmology and Astrophysics from Relaxed Galaxy Clusters III: Thermodynamic Profiles and Scaling Relations
 
To improve
----------

 * use M200c instead to be consistent with WL
 * weak lensing, M500
 * look at weighting the giant
 * M200 mean in cosmo codes
 * changing the slope as s function of mass
 * add substructure
 
"""
from scipy.stats import lognorm
from scipy.stats import norm
#import cPickle
#import fileinput
import astropy.io.fits as fits
import astropy.cosmology as co
import astropy.units as u
import astropy.constants as cc
cosmo = co.Planck13
from scipy.interpolate import interp1d
from os.path import join
import os
import astropy.units as uu
import numpy as n
import glob
import scipy.spatial.ckdtree as t
import time

class ClusterScalingRelations_Mantz2016() :
	"""
	Loads the environement to assign stellar masses to halos from dark matter only simulations, here MultiDark simulations.
	:param Lbox: length of the box in Mpc/h 
	:param wdir: Path to the multidark lightcone directory
	:param boxDir: box directory name
	:param snl: list of snapshots available
	:param zsl: list of redshift corresponding to the snapshots   
	:param zArray: redshift array to be considered to interpolate the redshift -- distance conversion
	:param Hbox: Hubble constant at redshift 0 of the box
	:param Melement: Mass of the resolution element in solar masses.   
	:param columnDict: dictionnary to convert column name into the index to find it in the snapshots
	
	definitions
	-----------
		- Planck flat LCDM cosmology
		- :math:`m = ln(M_{500} / (10^{15} M_\odot))`
		- :math:`m_{gas} = ln(M_{gas, 500} / (10^{15} M_\odot))` is the gas mass within r500
		- :math:`m_{lens} = ln(M_{lens, 500} / (10^{15} M_\odot))` is the spherical mass estimate from lensing corresponding to an
	idealized shear profile without statistical noise
		- :math:`l = ln(L_{500} / (E(z) 10^{44} erg/s))` where L is defined as the cluster rest-frame luminosity in the 0.1 - 2.4 keV band.
		- :math:`t = ln(kT_{500} / keV)` where kT is the emission weighted temperature measured in annulus from 0.15 to 1.0 r500
		- :math:`E(z) = H(z)/H_0`
		- :math:`\epsilon = ln(E(z))`
	
	Parameters
	----------
	
		* normalization for parameter X, :math:`N_X`
		* slope for E(z) for parameter X, :math:`slope_E_X`
		* slope for M500 for parameter X, :math:`slope_{M500}_X`
	 
	
	Workflow
	--------
		- Select relaxed clusters from the DM point of view. to be defined how ... with T/U ? Look at the publications from Sembolini, Yepes, Knebe ...
		- using M15, add gas density profile and temperature using scaling relations
	
	"""

	def __init__(self,Lbox=2500.0 * uu.Mpc, boxDir=os.environ['MD04'], snl=[], Hbox = 67.77 * uu.km / (uu.s * uu.Mpc), Melement = 23593750000.0 ):
		self.Lbox = Lbox # box length
		self.Hbox = Hbox # Hubble constant at redshift 0 in the box
		self.boxDir = boxDir # directory of the box where the snapshots a stored
		self.snl = snl # snapshot list, path to files
		self.Melement = Melement # mass of one particle in the box
		self.h = 0.6777
		# parameters used to run the simulation
		self.omega_lambda = 0.692885
		self.omega_matter = 0.307115
		self.omega_baryon = 0.048206
		self.ns = 0.96
		self.sigma8 = 0.8228
		self.G = 6.67428 * 10**(-9) # cm3 g-1 s-2
		self.Msun = 1.98892 * 10**(33.) # g
		self.Npart = 3840
		self.force_resolution = 5. # kpc /h
		
		self.N_Mgas = 31.98
		self.N_kT 	= 2.18
		self.N_L 	= 103.7
		self.N_Lce 	= 102.66
	
		self.slope_E_Mgas 	= -0.11
		self.slope_E_kT 	= 0.61
		self.slope_E_L 		= 1.20
		self.slope_E_Lce 	= 1.82

		self.slope_M500_Mgas= 1.04
		self.slope_M500_kT 	= 0.66
		self.slope_M500_L 	= 1.26
		self.slope_M500_Lce = 1.06
		
		self.E035 = cosmo.efunc(0.35)
		
		# converts logM500 to clusters observables
		self.m500_to_qty = lambda logM500, z, slope_efunc, slope_m500, normalization : n.e**normalization * (cosmo.efunc(z)/self.E035)**(slope_efunc) * (10**(logM500-n.log10(6)-14))**(slope_m500)
		
		self.logM500_to_logMgas = lambda logM500, z : self.m500_to_qty( logM500, z, self.slope_E_Mgas, self.slope_M500_Mgas, self.N_Mgas)
		self.logM500_to_kT 		= lambda logM500, z : self.m500_to_qty( logM500, z, self.slope_E_kT, self.slope_M500_kT, self.N_kT)
		self.logM500_to_L 		= lambda logM500, z : self.m500_to_qty( logM500, z, self.slope_E_L, self.slope_M500_L, self.N_L)
		self.logM500_to_Lce		= lambda logM500, z : self.m500_to_qty( logM500, z, self.slope_E_Lce, self.slope_M500_Lce, self.N_Lce)
	
class ClusterScalingRelations_Zandanel14() :
	"""
	Loads the environement to assign stellar masses to halos from dark matter only simulations, here MultiDark simulations.
	:param Lbox: length of the box in Mpc/h 
	:param wdir: Path to the multidark lightcone directory
	:param boxDir: box directory name
	:param snl: list of snapshots available
	:param zsl: list of redshift corresponding to the snapshots   
	:param zArray: redshift array to be considered to interpolate the redshift -- distance conversion
	:param Hbox: Hubble constant at redshift 0 of the box
	:param Melement: Mass of the resolution element in solar masses.   
	:param columnDict: dictionnary to convert column name into the index to find it in the snapshots
	"""

	def __init__(self,Lbox=2500.0 * uu.Mpc, boxDir=os.environ['MD04'], snl=[], Hbox = 67.77 * uu.km / (uu.s * uu.Mpc), Melement = 23593750000.0 ):
		self.Lbox = Lbox # box length
		self.Hbox = Hbox # Hubble constant at redshift 0 in the box
		self.boxDir = boxDir # directory of the box where the snapshots a stored
		self.snl = snl # snapshot list, path to files
		self.Melement = Melement # mass of one particle in the box
		self.h = 0.6777
		# parameters used to run the simulation
		self.omega_lambda = 0.692885
		self.omega_matter = 0.307115
		self.omega_baryon = 0.048206
		self.ns = 0.96
		self.sigma8 = 0.8228
		self.G = 6.67428 * 10**(-9) # cm3 g-1 s-2
		self.Msun = 1.98892 * 10**(33.) # g
		self.Npart = 3840
		self.force_resolution = 5. # kpc /h
	
	def logkbTci(self, Mh_500, z, h=0.7, h70=1.0):
		"""
		Computes the cluster temperature no centrally excised after Mantz et al. 2010b.

		Returns : log of the temperature in keV
		
		:param Mh_500: halo mass array
		:param z: redshift array
		"""
		return 0.91 + 0.46 * n.log10(cosmo.efunc(z) * MH_500 * h / (h70 * 10**15 ) )

	def get_Tci_K(self, Mh_500, z, h=0.7, h70=1.0):
		"""
		Samples with a Gaussian the Tci -- M500 scaling relation

		returns the temperature in Kelvin
		"""
		return 10**(norm.rvs(loc = self.logkbTci(Mh_500, z, h=0.7, h70=1.0), scale = 0.06 ) )/cc.k_B.to('keV/K').value

	def get_kbTci_keV(self, Mh_500, z, h=0.7, h70=1.0):
		"""
		Samples with a Gaussian the Tci -- M500 scaling relation

		returns the temperature in keV
		"""
		return 10**(norm.rvs(loc = self.logkbTci(Mh_500, z, h=0.7, h70=1.0), scale = 0.06 ) )

	def logLxbol(self, Mh_500, z, h70=1.0):
		"""
		Computes the bolometric X-ray luminosity for a halo of mass Mh_500
		Returns : the bolometrix X-ray luminosity in keV
		
		:param Mh_500: halo mass array
		:param z: redshift array
		"""
		return -21.5 + 1.5 * n.log10(cosmo.efunc(z) * MH_500 * self.h / (h70 * 10**15 ) )
