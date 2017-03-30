
"""
.. class:: MultiDark

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class MultiDark is a wrapper to handle Multidark simulations results / outputs.

"""
from scipy.stats import lognorm
from scipy.stats import norm
import cPickle
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
		aexp = z/(1+z)
		return 2 * ( 0.0351 - 0.0247 * aexp) / ((Mh/ (10**(11.59 + 1.195 * aexp)) )**(- 1.376 + 0.826 * aexp) + ( Mh /(10**(11.59 + 1.195 * aexp)) )**(0.608 + 0.329 *aexp) ) #- 0.0225

	def meanSM(self, Mh, z):
		"""
		Computes the mu star parameter for a halo mass according to the Moster et al. 2013 equations
		
		Returns :
		
		$\mu_* = 2\left(0.0351 - 0.0247 \frac{z}{1+z}\right)\left(\left[\left(\frac{10^{11.59 + 1.195 \frac{z}{1+z}}}{M_h}\right)^{1.376 - 0.826 \frac{z}{1+z}} + \left(\frac{M_h}{10^{11.59 + 1.195 \frac{z}{1+z}}} \right)^{0.608 + 0.329 \frac{z}{1+z}}  \right]^{-1}\right)- 0.0225$
		
		
		:param Mh: halo mass array
		:param z: redshift array
		"""
		aexp = z/(1+z)
		return n.log10(Mh * 2 * ( 0.0351 - 0.0247 * aexp) / ((Mh/ (10**(11.59 + 1.195 * aexp)) )**(- 1.376 + 0.826 * aexp) + ( Mh /(10**(11.59 + 1.195 * aexp)) )**(0.608 + 0.329 *aexp) )) #- 0.0225
	
	def sample_Ms( self, Mh, z, scatter = 0.15 ):
		"""
		Draws a stellar mass from a lognormal distribution centered on mu_star with witdth sigma_star
		
		:param Mh: halo mass
		:param z: redshift
		:param scatter: scatter in the stellar mass to halo mass relation
		"""
		return norm.rvs( loc = self.meanSM(Mh, z), scale = scatter )
		
	      
	      