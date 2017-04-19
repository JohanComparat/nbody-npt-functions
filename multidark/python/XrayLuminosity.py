
"""
.. class:: XrayLuminosity

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class XrayLuminosity is a wrapper to add Xray luminosity to the Multidark simulations results / outputs.

It adds Xray luminosity to simulations following the Bongiorno et al. 2016 model.

See http://adsabs.harvard.edu/abs/2016A%26A...588A..78B
equations 7, 10, 11, 12


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


class XrayLuminosity() :
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
		
		# parameters for the model
		self.z0 = 1.1
		self.psiStar = - 6.86
		
		
	def fz ( self, z ): 
		"""
		Computes the redshift component of the model $f_z(z)$ (equation 12)  
		
		:param z: redshift array
		:param z0: redshift turn over
		
		"""
		return n.piecewise(z, [z <= self.z0, z > self.z0], [ lambda z : (1.+z)**(5.82), lambda z : (1. + self.z0)**(5.82) * ((1.+z)/(1.+self.z0))**(2.36)])
	
        
	def fM(self, logM, z):
		"""
		Computes stellar mass component of the model $f_*$ (equation 10)
		
		:param logM: stellar mass array
		:param z: redshift array
		"""
		return (10**logM / 10**10.99 )**(0.24)*n.e**( - 10**logM / 10**10.99 )
	
	def fll( self, logM, z, ll ):
		"""
		Computes the specific accretion rate component of the model $f_\lambda$ (equation 11)
		:param logM: stellar mass
		:param z: redshift
		:param ll: log lambda SAR, specific accretion rate
		"""
		ll0 = 10**(33.8 - 0.48 * (logM - 11.) )
		ll_val = 10**ll
		g1z = 1.01 - 0.58 * (z - self.z0)
		#g2 = -3.72
		return ( ((ll_val)/(ll0))**(g1z) + ((ll_val)/(ll0))**(3.72) )**(-1.)
		
	def psi(self, logM, z, ll):
		"""
		Computes the bivariate distribution function (equation 7)
		:param logM: stellar mass
		:param z: redshift
		:param ll: log lambda SAR, specific accretion rate
		"""
		return 10**self.psiStar * self.fll( logM, z, ll ) * self.fM( logM, z ) * self.fz( z )
	      
