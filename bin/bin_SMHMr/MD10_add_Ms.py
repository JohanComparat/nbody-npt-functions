# overall python packages
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys 
# specific functions
from scipy.stats import norm
# dedicated packages
import StellarMass
sm = StellarMass.StellarMass()

def create_catalogs_out(fileList, z, minMS = 9.7):
	"""
	Adds stellar mass using the Moster et al. 2013 model to the rockstar outputs. 
	"""
	for fileName in fileList:
		t0=time.time()
		outFile = fileName[:-5]+"_Ms.fits"
		hd = fits.open(fileName)
		Mgal_mvir_Mo13 = norm.rvs( loc = sm.meanSM(10**hd[1].data['mvir'], z), scale = 0.15 )-n.log10(0.6777)
		sel = (Mgal_mvir_Mo13>minMS)
		
		col00 = fits.Column(name='stellar_mass_Mo13_mvir',format='D', unit='logMsun', array = Mgal_mvir_Mo13 )
		col01 = fits.Column(name='stellar_mass_reliable', format='L', array = sel )

		#define the table hdu 
		colArray = []
		colArray.append(hd[1].columns[0])
		# Mvir stellar mass
		colArray.append(col00)
		colArray.append(col01)

		hdu_cols  = fits.ColDefs(colArray)
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )

		#define the header
		prihdr = fits.Header()
		prihdr['author'] = 'JC'
		prihdr['SAMfile'] = os.path.basename(fileName)
		prihdr['minMS'] = minMS
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		if os.path.isfile(outFile):
			os.system("rm "+outFile)

		thdulist.writeto(outFile)
		print time.time()-t0

# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for ii in range(len(summ)):
	print summ[ii]
	fileList = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+summ['snap_name'][ii]+'_SAM_Nb_?.fits')))
	#outFile = fileName[:-5]+"_Ms.fits"
	z = summ['redshift'][ii]
	print fileList
	create_catalogs_out(fileList, z)

#create_catalogs(aexp = 0.75440, env='MD25', file_type="hlist")
sys.exit()


def create_catalogs(aexp = 0.74230, env='MD10' , file_type= "hlist", dV=-9999):
	"""
	old version
	"""
	
	if env=='MD04' :
		minMS = 7.2
	if env=='MD10' :
		minMS = 9.7
	if env=='MD25' :
		minMS = 11.3

	fileList = n.array(glob.glob(os.path.join(os.environ[env], "snapshots", file_type+"_*_SAM_Nb_*.fits" )))
	fileList.sort()
	z = 1./0.74230 -1.
	fileList.sort()
	print fileList
	for fileName in fileList:
		t0=time.time()
		outFile = os.path.join(os.environ[env], "catalogs", os.path.basename(fileName)[:-5] + "_stellar_mass.fits")
		print outFile
		hd = fits.open(fileName, mode='update')
		
		Mgal_mvir_Mo13 = norm.rvs( loc = sm.meanSM(10**hd[1].data['mvir'], z), scale = 0.15 )-n.log10(0.6777)
		sel = (Mgal_mvir_Mo13>minMS)
		
		hd.flush() 
		hd.close()
		VAC	
		col00 = fits.Column(name='stellar_mass_Mo13_mvir',format='D', unit='logMsun', array = Mgal_mvir_Mo13 )
		col01 = fits.Column(name='stellar_mass_reliable', format='L', array = sel )
		
		#define the table hdu 
		colArray = []
		for col in hd[1].columns :
			colArray.append(col)
		
		# Mvir stellar mass
		colArray.append(col00)
		colArray.append(col01)
		
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

hde = fits.open('plates-dr14.fits')
hd2 = fits.open('test.fits')

hdr = fits.Header()
hdr.append(('stellar_mass_Mo13_mvir', 'float'))
hdr.append(('stellar_mass_reliable', 'boolean'))

Mgal_mvir_Mo13 = n.ones_like(hd[1].data['PLATE'])
data = n.transpose([ Mgal_mvir_Mo13, Mgal_mvir_Mo13+1] )

# fits update 
fits.update('test.fits', data, hdr, 2)


colArray = []
for col in hd[1].columns :
	colArray.append(col)

# Mvir stellar mass
colArray.append(col00)

hdu_cols  = fits.ColDefs(colArray)

tb_hdu = fits.Binhdu	TableHDU.from_columns( hdu_cols )

# fits update 


fits.append('test.fits', n.array([n.ones_like(hd[1].data['PLATE']), n.ones_like(hd[1].data['PLATE'])*2.], dtype=(n.record, [('col1', '>f8'), ('col2', '>f8')])), hdr, 2)

hd3 = fits.open('test.fits')
hd3