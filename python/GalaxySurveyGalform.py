"""
.. class:: GalaxySurveyGalform

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class GalaxySurveyGalform is dedicated to handling Galform fake galaxy mock catalog.

"""
from os.path import join
import os
import numpy as n
import astropy.io.fits as fits
from MiscellanousFunctionsLibrary import *
import astropy.units as u

class GalaxySurveyGalform:
	"""
	Loads the environement proper to the Galform mock catalog:
	 * Defines all the proper paths in the database,
         * Opens the catalog,
        :param redshift_catalog: name of the galform redshift catalog (path to the fits file)
	"""
	def __init__(self,redshift_catalog="galform.elg.fits"):
		self.redshift_catalog = redshift_catalog
		self.database_dir = os.environ['DATA_DIR']
		self.galform_dir = join(self.database_dir, "Simulations","galform-lightcone")
		self.galform_catalog_dir = join(self.galform_dir,"catalogs")
		hd = fits.open(join(self.galform_catalog_dir,self.redshift_catalog))
		self.catalog = hd[1].data
		hd.close()

	def computeLineLuminosity(self,line,distanceCorrection):
		""" computes the line luminosities for the line list given.
		:param catalog: fits catalog containing redshift, EBV and line fluxes
		:param line:
		"""
		#ebvCorrection=n.array([ 10**(0.4 *self.catalog['EBV_MW'][i] * CalzettiLaw((1 + self.catalog['zspec'][i]) * line[1])) for i in range(len(self.catalog['zspec']))])
		#correctionAperture = 1. / self.catalog['fo']
		#flux=ebvCorrection * correctionAperture * self.catalog[line[2]+'_flux']* u.erg/ u.cm**2 /u.s
		flux= self.catalog[line[2]+'_flux']* u.erg/ u.cm**2 /u.s
		Luminosity=fits.Column(name=line[2]+"_luminosity",format="D", unit="erg/s", array=distanceCorrection*flux )
		return Luminosity


