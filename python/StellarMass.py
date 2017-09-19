
"""
.. class:: MultiDark

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class MultiDark is a wrapper to handle Multidark simulations results / outputs.

"""
from scipy.stats import lognorm
from scipy.stats import norm
#import cPickle
import fileinput
import astropy.io.fits as fits
import astropy.cosmology as co
import astropy.units as u
c2 = co.Planck13
from scipy.interpolate import interp1d
from os.path import join
import os
import astropy.units as uu
import numpy as n
import glob
import scipy.spatial.ckdtree as t
import time

class StellarMass() :
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
	
	def SMHMr(self, Mh, z):
		"""
		Computes the mu star parameter for a halo mass according to the Moster et al. 2013 equations
		
		Returns :
		
		$\mu_* = 2\left(0.0351 - 0.0247 \frac{z}{1+z}\right)\left(\left[\left(\frac{10^{11.59 + 1.195 \frac{z}{1+z}}}{M_h}\right)^{1.376 - 0.826 \frac{z}{1+z}} + \left(\frac{M_h}{10^{11.59 + 1.195 \frac{z}{1+z}}} \right)^{0.608 + 0.329 \frac{z}{1+z}}  \right]^{-1}\right)- 0.0225$
		
		
		:param Mh: halo mass array
		:param z: redshift array
		"""
		aexp = z/(1.+z)
		return 2. * ( 0.0351 - 0.0247 * aexp) / ((Mh/ (10**(11.59 + 1.195 * aexp)) )**(- 1.376 + 0.826 * aexp) + ( Mh /(10**(11.59 + 1.195 * aexp)) )**(0.608 + 0.329 *aexp) ) #- 0.0225

	def meanSM(self, Mh, z):
		"""
		Computes the mu star parameter for a halo mass according to the Moster et al. 2013 equations
		
		Returns :
		
		$\mu_* = 2\left(0.0351 - 0.0247 \frac{z}{1+z}\right)\left(\left[\left(\frac{10^{11.59 + 1.195 \frac{z}{1+z}}}{M_h}\right)^{1.376 - 0.826 \frac{z}{1+z}} + \left(\frac{M_h}{10^{11.59 + 1.195 \frac{z}{1+z}}} \right)^{0.608 + 0.329 \frac{z}{1+z}}  \right]^{-1}\right)- 0.0225$
		
		
		:param Mh: halo mass array
		:param z: redshift array
		"""
		aexp = z/(1.+z)
		return n.log10(Mh * 2. * ( 0.0351 - 0.0247 * aexp) / ((Mh/ (10**(11.59 + 1.195 * aexp)) )**(- 1.376 + 0.826 * aexp) + ( Mh /(10**(11.59 + 1.195 * aexp)) )**(0.608 + 0.329 *aexp) )) #- 0.0225
	
	def sample_Ms( self, Mh, z, scatter = 0.15 ):
		"""
		Draws a stellar mass from a lognormal distribution centered on mu_star with witdth sigma_star
		
		:param Mh: halo mass
		:param z: redshift
		:param scatter: scatter in the stellar mass to halo mass relation
		"""
		return norm.rvs( loc = self.meanSM(Mh, z), scale = scatter )
		
class StellarMass_EMERGE_Moster2017() :
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

	def __init__(self,Lbox=1000.0 * uu.Mpc, boxDir=os.environ['MD10'], snl=[], Hbox = 67.77 * uu.km / (uu.s * uu.Mpc), Melement = 23593750000.0 ):
		# parameters related to the simulations
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
		# parameters and equations related to EMERGE
		
		# equation (7)
		self.log_M0 = 11.339 # +0.005 -0.080
		self.log_Mz = 0.692 # +0.010 -0.009
		self.log10_M1 = lambda z : self.log_M0  + self.log_Mz * (z/(1.+z)) 
		
		# equation (8)
		self.epsilon_0 = 0.005 
		self.epsilon_z = 0.689
		self.epsilon_N = lambda z : self.epsilon_0  + self.epsilon_z * (z/(1.+z)) 
		
		# equation (9)
		self.beta_0 = 3.334 
		self.beta_z = -2.079
		self.beta = lambda z : self.beta_0  + self.beta_z * (z/(1.+z)) 
		
		# equation (10)
		self.gamma_0 = 0.966
		self.gamma = lambda z : self.gamma_0
		
		# equation (5) <= (7, 8, 9, 10)
		# integrated efficiency function of mass and redshift
		self.epsilon = lambda stellar_mass, z : 2. * self.epsilon_N(z) /((stellar_mass / 10**self.log10_M1(z))**(-self.beta(z)) + (stellar_mass / 10**self.log10_M1(z))**(self.gamma(z)))
		
		# equation (6)
		# mass at which baryon conversion is most efficient
		self.M_max = lambda z : 10**self.log10_M1(z) * (self.beta(z)/self.gamma(z))**(1/(self.beta(z) + self.gamma(z)))
		
		# equation (13)
		self.tau_0 = 4.282 
		self.tau_s = 0.363
		self.tau = lambda t_dyn, stellar_mass : t_dyn * self.tau_0 * (stellar_mass * 10**(-10.))**(-self.tau_s)
		
		# equation (14), stripping
		self.f_s = 0.122
		# equation (15), merging
		self.f_esc = 0.338 
		
	def reconsitute_history(self):
		"""
		reads a fits file at a given redshift:
		#. split central - sat
		 #. read and match to its predecessors at the previous redshift for centrals. 
		 #. read and match at all previous redshifts for sat
		 #. 2 catalogs of matched properties
		#. write history catalogs with properties of interest
		
		#. retrieve the properties of interest
		#. 
		columns available in short files 
		'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'pid': 33
		"""
		return 0.
	
	def sample_stellar_mass(self):
		"""
		Given a file written by reconstitute history, 
		#. computes the galaxy properties 
		#. writes them to a new file "_galaxy.fits"
		"""
		return 0.
