
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

from scipy.integrate import quad

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
		
		# parameters of the model
		self.z0 = 1.1
		self.psiStar = - 6.86
		
		
	def f_z ( self, z ): 
		"""
		Computes the redshift component of the model :math:`f_z(z)` (equation 12)  
		
		:param z: redshift array
		:param z0: redshift turn over
		
		"""
		return n.piecewise(z, [z <= self.z0, z > self.z0], [ lambda z : (1.+z)**(5.82), lambda z : (1. + self.z0)**(5.82) * ((1.+z)/(1.+self.z0))**(2.36)])
	
        
	def f_Mstar(self, logM, z):
		"""
		Computes stellar mass component of the model $f_*$ (equation 10)
		
		:param logM: log stellar mass array
		:param z: redshift array
		"""
		return (10**(logM - 10.99) )**(0.24) * n.e**( - 10**(logM - 10.99) )
	
	def f_lambda_sar( self, logM, z, log_lambda_SAR ):
		"""
		Computes the specific accretion rate component of the model $f_\lambda$ (equation 11)
		:param logM: stellar mass
		:param z: redshift
		:param log_lambda_SAR: log lambda SAR, specific accretion rate
		"""
		log_lambda_SAR_var = 10**( log_lambda_SAR - 33.8 + 0.48 * (logM - 11.) )		
		g1z = 1.01 - 0.58 * (z - self.z0)
		return 1. / ( log_lambda_SAR_var**(g1z) + log_lambda_SAR_var**(3.72) )
		
	def psi_log(self, log_lambda_SAR, logM, z):
		"""
		Computes the bivariate distribution function (equation 7)
		:param logM: stellar mass
		:param z: redshift
		:param log_lambda_SAR: log lambda SAR, specific accretion rate
		"""
		return 10**self.psiStar * self.f_lambda_sar( logM, z, log_lambda_SAR ) * self.f_Mstar( logM, z ) * self.f_z( z )
	      
	def psi(self, lambda_SAR, stellar_mass, redshift):
		"""
		Computes the bivariate distribution function (equation 7)
		:param logM: stellar mass
		:param z: redshift
		:param log_lambda_SAR: log lambda SAR, specific accretion rate
		"""
		return psi_log(n.log10(lambda_SAR), n.log10(stellar_mass), redshift)
	
	def Phi_lambda_SAR(self, log_lambda_SAR, redshift):
		"""
		Integrates the bivariate distribution function on logM (equation 9) gives the specific accretion rate distribution function (SARDF). :math:`dN/dVdlogM`
		:param logM: stellar mass
		:param z: redshift
		:param log_lambda_SAR: log lambda SAR, specific accretion rate
		"""
		integrand = lambda M, log_lambda_SAR, z : self.psi_log( log_lambda_SAR, n.log10(M), z)/(n.log(10)*M)
		return quad(integrand, 10**9.5, 10**10.5, args=(log_lambda_SAR, redshift))[0] + quad(integrand, 10**10.5, 10**11.5, args=(log_lambda_SAR, redshift))[0] + quad(integrand, 10**11.5, 10**12.5, args=(log_lambda_SAR, redshift))[0]+quad(integrand, 10**12.5, 10**13.5, args=(log_lambda_SAR, redshift))[0]
	
	def Phi_stellar_mass(self, stellar_mass, redshift):
		"""
		Integrates the bivariate distribution function on logM (equation 8). Gives the host galxy stellar mass function (AGN HGMF). :math:`dN/dVdlog\lambda`
		:param logM: stellar mass
		:param z: redshift
		"""
		integrand = lambda lambda_SAR, logM, z : self.psi_log( n.log10(lambda_SAR), logM, z)/(n.log(10)*lambda_SAR)
		return quad(integrand, 10**32, 10**33, args=(stellar_mass, redshift))[0] + quad(integrand, 10**33, 10**34, args=(stellar_mass, redshift))[0] + quad(integrand, 10**34, 10**35, args=(stellar_mass, redshift))[0]+quad(integrand, 10**35, 10**36, args=(stellar_mass, redshift))[0]
	
	def Phi_stellar_mass_to_X(self, X, stellar_mass, redshift):
		"""
		Integrates the bivariate distribution function on logM (equation 8). Gives the host galxy stellar mass function (AGN HGMF). :math:`dN/dVdlog\lambda`
		:param logM: stellar mass
		:param z: redshift
		"""
		integrand = lambda lambda_SAR, logM, z : self.psi_log( n.log10(lambda_SAR), logM, z)/(n.log(10)*lambda_SAR)
		if X>35 :
			return quad(integrand, 10**32, 10**33, args=(stellar_mass, redshift))[0] + quad(integrand, 10**33, 10**34, args=(stellar_mass, redshift))[0] + quad(integrand, 10**34, 10**35, args=(stellar_mass, redshift))[0]+quad(integrand, 10**35, 10**35, args=(stellar_mass, redshift))[0]
		elif X>34 :
			return quad(integrand, 10**32, 10**33, args=(stellar_mass, redshift))[0] + quad(integrand, 10**33, 10**34, args=(stellar_mass, redshift))[0] + quad(integrand, 10**34, 10**X, args=(stellar_mass, redshift))[0]
		elif X>33 :
			return quad(integrand, 10**32, 10**33, args=(stellar_mass, redshift))[0] + quad(integrand, 10**33, 10**X, args=(stellar_mass, redshift))[0]
		else :
			return quad(integrand, 10**32, 10**X, args=(stellar_mass, redshift))[0]
	
	def obscured_fraction_optical_Merloni2015(self, logLX):
		"""
		Observed obscured fraction as a function of luminosity shown in Fig. 6 of Merloni et al. 2014 'the incidence of obscuration in active galactic nuclei'.
		Equation 1 wih the parameters given in page 3556
		"""
		return 0.56 + n.arctan((43.89 - logLX)/0.46)/n.pi
	
	def obscured_fraction_HRz_Merloni2015(self, logLX):
		"""
		Observed obscured fraction as a function of luminosity shown in Fig. 7 of Merloni et al. 2014 'the incidence of obscuration in active galactic nuclei'.
		
		No evolution with X ray luminosity, data consistent with a constant fraction.
		"""
		return 0.56
	