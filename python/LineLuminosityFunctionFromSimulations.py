"""
.. class:: LineLuminosityFunctionFromSimulations

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class LineLuminosityFunctionFromSimulations is dedicated to measuring the line luminosity functions obtained from simulations.
"""
from os.path import join
import os
import astropy.cosmology as co
cosmo=co.FlatLambdaCDM(H0=70,Om0=0.3)
import astropy.io.fits as fits
import numpy as n
from scipy.optimize import curve_fit


class LineLuminosityFunctionFromSimulations:
	"""
	The line luminosity function class
	:param lineWavelength: restframe wavelength in the air
	:param lineName: name of the line used in the catalogs.
	:param cosmology: cosmology used (astropy class) Default H0=70,Omega matter=0.3
    :param surveyName: Name of the survey used (needs to be the one given in the database)
    :param redshift_catalog: name of the redshift catalog 
	:param luminosityBins: bins in luminosity equally spaced in log space.
	:param outputFolder: folder where the results will be written
	:param zmin: minimum redshift included
	:param zmax: maximum redshift included
	"""
	def __init__(self, lineWavelength=3727.4228417998916, lineName="OII3727", cosmology = cosmo, surveyName ="GALFORM", surveyDir =  join("Simulations","galform-lightcone"), redshift_catalog = "galform.ELG.fits", luminosityBins = n.logspace(38,45,50), outputFolder="emissionLineLuminosityFunctions" , zmin=0.6, zmax=0.8):
		self.lineWavelength = lineWavelength
		self.lineName = lineName
		self.cosmology = cosmology
		self.surveyName = surveyName
		self.redshift_catalog = redshift_catalog
		self.database_dir = os.environ['DATA_DIR']
		self.survey_dir = join(self.database_dir , surveyDir)
		self.catalog_dir = join(self.survey_dir,"catalogs")
		self.output_dir = join(self.survey_dir,"products",outputFolder,lineName)
		os.system('mkdir '+self.output_dir)

		hd = fits.open(join(self.catalog_dir,self.redshift_catalog))
		self.catalog = hd[1].data
		hd.close()

		self.Ngalaxies = len(self.catalog)
		#self.nbins = 15#n.arange(38.5,45,0.25)#15
		self.luminosityBins = luminosityBins #15
		#self.nbinsUD = 4

		self.zmin = zmin
		self.zmax = zmax

		self.luminosity = self.catalog[lineName+'_luminosity']
		self.volume_per_sq_degree=lambda z1,z2 : (cosmo.comoving_volume( z2 ) - cosmo.comoving_volume( z1 )) *n.pi/129600.


	def setRedshiftArray(self,redshiftColumn='zObs'):
		""" sets the redshift array
		:param redshiftColumn: column of the catalog corresponding to the redshift.
		Stores it in self.redshift.
		"""
		self.redshift = self.catalog[redshiftColumn]

	def setRedshiftSelection(self):
		""" sets the redshift selection
		:param redshiftQualityColumn: column of the catalog corresponding to the  quality of the redshifts.
		:param lowerBound : lower bound to redshift quality :  zquality > lowerBound
		:param upperBound : upper bound to the redshift quality : zquality < upperBound
		Stores it in self.redshiftSelection.
		"""
		self.redshiftSelection = ( self.redshift>self.zmin ) & ( self.redshift<self.zmax )

	def setWeightArray(self,weightColumn):
		""" sets the weight column 
		:param weightColumn: statistical weight per galaxy 1 / (area * TSR * SSR)
		Divides the weight by the volume of the bin stores it in self.weight.
		"""
		self.weight = n.ones_like(self.luminosity) * weightColumn / self.volume_per_sq_degree(self.zmin,self.zmax)

	def computeMeanWeightedRedshift(self,sel):
		""" Computes the weighted mean redshift of the sample.
		"""
		selection = (sel) & (self.redshiftSelection)
		self.meanRedshift = n.average(self.redshift[selection], weights = self.weight[selection])

	def computeHistogramLF(self,sel):
		""" Computes the weighted and unweighted histogram to get the number density and Poisson errors.
		:param sel: array selecting the galaxies of interest in the catalog (Boolean). 
		Returns Weighted density, Error on the weighted density, Number of galaxies used in eacah bin, the luminosity bins.
		It stores the values in self.LF, self.LFerr_poisson, self.ngals. It also evaluates the mean luminosity in each luminosity bin self.xL and dlogL to obtain the LF
		"""
		selection = (sel) & (self.redshiftSelection)
		N10p,bin1p=n.histogram(self.luminosity[selection],bins=self.luminosityBins)
		N10,bin1=n.histogram(self.luminosity[selection], bins= self.luminosityBins, weights= self.weight[selection] )
		self.LF, self.LFerr_poisson, self.ngals = N10, N10*N10p**0.5/N10p, N10p

		xSelections=n.array([ (self.luminosity > self.luminosityBins[ii]) &(self.luminosity< self.luminosityBins[ii+1] ) & (selection) for ii in range( len( self.luminosityBins ) -1 ) ])

		xLi= []
		for jj in range(len(xSelections)) :
			if len(self.luminosity[xSelections[jj]])>0:
				xLi.append( n.average( self.luminosity[xSelections[jj]], weights= self.weight[xSelections[jj]] ) )
			else:
				xLi.append( (self.luminosityBins[jj]+self.luminosityBins[jj+1])/2. )
		
		self.xL=n.array(xLi)
		dLogL_all = (self.luminosityBins[1:] - self.luminosityBins[:-1]) / ((self.luminosityBins[1:] + self.luminosityBins[:-1])/2.)
		self.dLogL = dLogL_all[0]

	def computeHistogramVariance(self,sel,jk=0.1):
		""" Computes the variance of the histogram using N subsamples.
		:param sel: array selecting the galaxies of interest in the catalog (Boolean). 
		:param jk: percentage of the data set removed in each realization.
		Stores the values in self.LFerr_jackknife
		"""
		selection = (sel) & (self.redshiftSelection)
		#N10p,bin1p=n.histogram(self.luminosity[selection],bins=self.luminosityBins)
		L_jk = self.luminosity[selection]
		w_jk = self.weight[selection]
		rdArr=n.random.rand(len(L_jk))
		values=n.arange(0,1+0.9*jk,jk)
		randSelNot=n.array([(rdArr>values[jj])&(rdArr<values[jj+1]) for jj in range(len(values)-1)])
		randSel=n.array([(el==False) for el in randSelNot])
		lumJK=[]
		for selR in randSel :
			N10,bin1=n.histogram(L_jk[selR], bins= self.luminosityBins, weights= w_jk[selR] )
			lumJK.append(N10)

		self.LFerr_jackknife = n.std(lumJK,axis=0)

	def get_completness_limit(self,sel):
		selection = (sel) & (self.redshiftSelection)
		bins=n.logspace(1,3,20)
		aa,bb = n.histogram(self.catalog[self.lineName+'_EW'][selection], bins=bins)
		self.completness_limit_EW = bb[n.argmax(aa)+3]
		
		EWselection = (self.catalog[self.lineName+'_EW'][selection] >0.9* self.completness_limit_EW )&( self.catalog[self.lineName+'_EW'][selection]<1.1* self.completness_limit_EW)

		self.completness_limit_luminosity = n.median( self.catalog[ self.lineName+'_luminosity'][ selection ][ EWselection ])

		# bins=n.logspace(39.5,43,20)
		# aa,bb = n.histogram(self.catalog[self.lineName+'_luminosity'][selection], bins=bins)
		#self.completness_limit_luminosity = bb[n.argmax(aa)+1]
	

	def writeLF(self,sel,surveyNameSuffix=""):
		""" writes the measured LF and the data used to derive it to an ascii and a fits file.
		"""
		filename = self.lineName + "-" + self.surveyName+surveyNameSuffix + "-z" + str( n.round( self.meanRedshift ,3 )) 

		selection = (sel) & (self.redshiftSelection)
		new_columns = self.catalog.columns
		hdu2 = fits.BinTableHDU.from_columns(new_columns)
		hdu2.data = hdu2.data[selection]
		hdu2.header.add_comment(str(self.completness_limit_luminosity))
		os.system('rm -rf '+ join(self.output_dir , filename + ".fits"))
		hdu2.writeto(join(self.output_dir , filename + ".fits"))                        
                        
		head= " Lmin Lmax Lmean phi phiErr_jk phiErr_poisson Ngalaxy"
		f=open(join(self.output_dir , filename + ".txt"),'w')
		n.savetxt(f, n.transpose([self.luminosityBins[:-1], self.luminosityBins[1:], self.xL, self.LF/self.dLogL, self.LFerr_poisson/self.dLogL, self.LFerr_jackknife /self.dLogL, self.ngals]) ,header= head)
		f.close()



