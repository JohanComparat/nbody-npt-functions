
"""
.. class:: ClusterScalingRelations

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class clusterScalingRelations is a wrapper to add cluster physics to the Multidark simulations results / outputs.

Based on 
 * Mantz et al. 2010b
 * Zandanel et al. 2014
"""
from scipy.stats import lognorm
from scipy.stats import norm
import cPickle
import fileinput
import astropy.io.fits as fits
import astropy.cosmology as co
import astropy.units as u
import astropy.constants as cc
c2 = co.Planck13
from scipy.interpolate import interp1d
from os.path import join
import os
import astropy.units as uu
import numpy as n
import glob
import scipy.spatial.ckdtree as t
import time

class ClusterScalingRelations() :
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
		return 0.91 + 0.46 * n.log10(c2.efunc(z) * MH_500 * h / (h70 * 10**15 ) )

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
		return -21.5 + 1.5 * n.log10(c2.efunc(z) * MH_500 * self.h / (h70 * 10**15 ) )
