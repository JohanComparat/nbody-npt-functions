
"""
.. class:: DarkSkies

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class DarkSkies is a wrapper to handle DarkSkies simulations results / outputs.

"""
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

class DarkSkiesSimulation :
	"""
	Loads the environement proper to the DarkSkies simulations. This is the fixed framework of the simulation.
			
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

	def __init__(self,Lbox=8000.0 * uu.Mpc, wdir = join(os.environ['DATA_DIR'],"DarkSkies"), boxDir="snapshots", snl=n.array(glob.glob(join(os.environ['DATA_DIR'],"DarkSkies", "snapshots","ds14_catalog_300particles.dat"))), zsl=None, zArray=n.arange(0.2,2.4,1e-1), Hbox = 67.77 * uu.km / (uu.s * uu.Mpc), Melement = 3.9*10**(10.0) ):
		self.Lbox = Lbox # box length
		self.Hbox = Hbox # Hubble constant at redshift 0 in the box
		self.wdir = wdir # working directory
		self.boxDir = boxDir # directory of the box where the snapshots a stored
		self.snl = snl # snapshot list
		self.zsl = zsl # corresponding redshift list
		self.zArray = zArray # redshift for the dC - z conversion
		self.Melement = Melement # mass of one particle in the box
		self.h = 0.6881
		self.ns = 0.9676
		self.G = 6.67428 * 10**(-9) # cm3 g-1 s-2
		self.Msun = 1.98892 * 10**(33.) # g
		self.Npart = 10240
		self.force_resolution = 36.8 # kpc /h
		self.columnDict = {'mvir': 0, 'vmax': 1, 'x': 2, 'y': 3, 'z': 4, 'id': 5, 'pid': 6}
		self.Melement = 3.9*10**(10.0) 
			
	def writePositionCatalogPM(self, ii, vmin=30., mmin=300*3.9*10**(10.0) , NperBatch = 10000000):
		"""
		Extracts the positions and velocity out of a snapshot of the Multidark simulation.        
		:param ii: index of the snapshot in the list self.snl
		:param vmin: name of the quantity of interest, mass, velocity.
		:param vmax: of the quantity of interest in the snapshots.
		:param NperBatch: number of line per fits file, default: 1000000
		 """		
		fl = fileinput.input(self.snl[ii])
		nameSnapshot = os.path.basename(self.snl[ii])[:-4]
		Nb = 0
		count = 0
		output = n.zeros((NperBatch,7))
		for line in fl:
			if line[1] == "#" :
				continue

			line = line.split()
			newline =n.array([int(line[self.columnDict['id']]), float(line[self.columnDict['pid']]), float(line[self.columnDict['x']]), float(line[self.columnDict['y']]), float(line[self.columnDict['z']]), float(line[self.columnDict['vmax']]), n.log10(float(line[self.columnDict['mvir']])) ])
			if   float(line[self.columnDict['vmax']])>vmin and float(line[self.columnDict['mvir']])>mmin :
				output[count] = newline
				count+=1
				
			if count == NperBatch  :
				#print "count",count
				#print output
				#print output.shape
				#print output.T[0].shape
				#define the columns
				col0 = fits.Column(name='id',format='D', array= output.T[0] )
				col1 = fits.Column(name='pid',format='D', array= output.T[1] )
				col2 = fits.Column(name='x',format='D', array=output.T[2] )
				col3 = fits.Column(name='y',format='D', array= output.T[3] )
				col4 = fits.Column(name='z',format='D', array= output.T[4] )
				col5 = fits.Column(name='vmax',format='D', array= output.T[5] )
				col6 = fits.Column(name='mvir',format='D', array=output.T[6] )
				#define the table hdu 
				hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				#define the header
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = nameSnapshot
				prihdr['count'] = count
				prihdr['batchN'] = Nb
				prihdr['author'] = 'JC'
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				os.system("rm "+self.snl[ii][:-4]+"_PM_Nb_"+str(Nb)+".fits")
				thdulist.writeto(self.snl[ii][:-4]+"_PM_Nb_"+str(Nb)+".fits")
				Nb+=1
				count=0
				#resest the output matrix
				output = n.zeros((NperBatch,7))
		
		
		# and for the last batch :		
		col0 = fits.Column(name='id',format='D', array= output.T[0][:count] )
		col1 = fits.Column(name='pid',format='D', array= output.T[1][:count] )
		col2 = fits.Column(name='x',format='D', array=output.T[2][:count] )
		col3 = fits.Column(name='y',format='D', array= output.T[3][:count] )
		col4 = fits.Column(name='z',format='D', array= output.T[4][:count] )
		col5 = fits.Column(name='vmax',format='D', array= output.T[5][:count] )
		col6 = fits.Column(name='mvir',format='D', array=output.T[6][:count] )
		#define the table hdu 
		hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		#define the header
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = nameSnapshot
		prihdr['count'] = count
		prihdr['batchN'] = Nb
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		os.system("rm "+self.snl[ii][:-4]+"_PM_Nb_"+str(Nb)+".fits")
		thdulist.writeto(self.snl[ii][:-4]+"_PM_Nb_"+str(Nb)+".fits")

	def computeSingleDistributionFunctionJKresampling(self, fileList, rootname, name, bins, Ljk = 100., overlap = 1. ) :
		"""
		Extracts the distribution of quantity 'name' out of all snapshots of the Multidark simulation.
		Resamples the box in smaller boxes of length Ljk in Mpc/h
		:param ii: index of the snapshot in the list self.snl
		:param name: name of the quantity of interest, mass, velocity.
		:param index: of the quantity of interest in the snapshots.
		:param bins: binning scheme to compute the historgram.
		:param Ljk: length of the resampled box
		:param overlap: allowed overlap between resampled realizations : 1 = no overlap 2 : 50% overlap ... 
		"""		
		output_dir = join(self.wdir,"properties",name)
		os.system('mkdir '+ output_dir)
		# define boundaries
		NBoundariesPerSide = int(overlap*self.Lbox.value/Ljk)
		bounds = n.arange(NBoundariesPerSide+1)* Ljk / overlap
		#print "boundaries on each side: ", bounds
		Xi, Yi, Zi = n.meshgrid(bounds[:-1],bounds[:-1],bounds[:-1])
		X = n.ravel(Xi)
		Y = n.ravel(Yi)
		Z = n.ravel(Zi)	
		#print X.min(), X.max(), len(X),len(bounds)
		# loops over the fileList : fits files with the data
		nnC = n.zeros((len(fileList),len(X),len(bins)-1))
		nnS = n.zeros((len(fileList),len(X),len(bins)-1))
		for jj, file in enumerate(fileList):
			#print file
			dd = fits.open(file)[1].data
			cen = (dd['pid']==-1)
			sat = (cen==False) # (dd['pid']>=1)
			#computes the histogram for each resampling of the file
			for ii, xel in enumerate(X):
				#print ii
				xmin, ymin, zmin, xmax, ymax, zmax = X[ii], Y[ii], Z[ii], X[ii]+Ljk, Y[ii]+Ljk, Z[ii]+Ljk
				sel = (dd['x']>=xmin)&(dd['x']<xmax)&(dd['y']>=ymin)&(dd['y']<ymax)&(dd['z']>=zmin)&(dd['z']<zmax)&(dd[name]>bins[0])&(dd[name]<bins[-1])
				#print len(dd[name][(sel)&(cen)]), len(dd[name][(sel)&(sat)])
				if len(dd[name][(sel)&(cen)])>=1:
					nnC[jj][ii] = n.histogram(dd[name][(sel)&(cen)], bins = bins)[0]
				if  len(dd[name][(sel)&(sat)])>=1:
					nnS[jj][ii] = n.histogram(dd[name][(sel)&(sat)], bins = bins)[0]
			
		f = open(join(output_dir, rootname +"_Central_JKresampling.pkl"),'w')
		cPickle.dump(n.sum(nnC,axis=0),f)
		f.close()
		f = open(join(output_dir,rootname +"_Satellite_JKresampling.pkl"),'w')
		cPickle.dump(n.sum(nnS,axis=0),f)
		f.close()
		n.savetxt(join(output_dir,rootname+"_"+name+"_JKresampling.bins"),n.transpose([bins]))

