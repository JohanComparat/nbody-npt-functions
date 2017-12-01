#import StellarMass
import XrayLuminosity

import numpy as n
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

meanSM= lambda Mh, z : n.log10(Mh * 2. * ( 0.0351 - 0.0247 * z/(1.+z)) / ((Mh/ (10**(11.79 + 1.5 * z/(1.+z))) )**(- 0.9 + 0.5  * z/(1.+z)) + ( Mh /(10**(11.79 + 1.5 * z/(1.+z))) )**(0.67 + 0.2 * z/(1.+z)) ) )

fun = lambda mmm : norm.rvs( loc = mmm, scale = 0.15 )

plotDir = os.path.join(os.environ['HOME'], 'wwwDir', "eRoMok", "snapshot_sanity_checks")
mhalos = n.arange(10,15,0.1)

def plot_SMHMR(sf, env='MD10'):
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", "out_"+el["snap_name"]+"_SAM_Nb_?_Ms.fits")))
	fileList.sort()
	snapList = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", "out_"+el["snap_name"]+"_SAM_Nb_?.fits")))
	snapList.sort()
	print fileList, snapList
	for fileMs, fileN in zip(fileList, snapList):
		print fileN, fileMs
		hM = fits.open(fileMs)
		hN = fits.open(fileN)
		mass = hM[1].data['stellar_mass_Mo13_mvir']
		mvir = hN[1].data['mvir'] - n.log10(0.6777)
		smhmr = mass - mvir
		selection = (hM[1].data['stellar_mass_reliable'])
		
		p.figure(1, (6,6))
		p.plot(mvir[selection], smhmr[selection], 'k,', rasterized = True, alpha = 0.5 )
		p.plot(mhalos, meanSM(10**mhalos, sf['redshift'])-mhalos , 'r', rasterized = True )
		p.plot(mhalos, meanSM(10**mhalos, sf['redshift'])-mhalos - 0.15, 'r--', rasterized = True )
		p.plot(mhalos, meanSM(10**mhalos, sf['redshift'])-mhalos + 0.15, 'r--', rasterized = True )
		p.xlabel('mvir')
		p.ylabel('stellar mass - mvir')
		#p.yscale('log')
		p.title(str(el['redshift']))
		p.grid()
		p.savefig(os.path.join(plotDir, os.path.basename(fileMs)[:-5] + "_smhmr_moster2013.jpg"))
		p.clf()


	
# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for el in summ:
	print el
	plot_SMHMR(el, env='MD10')
