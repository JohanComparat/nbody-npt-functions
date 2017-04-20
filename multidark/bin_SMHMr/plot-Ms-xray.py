import StellarMass
import XrayLuminosity

import numpy as n
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

print " set up box, and redshift "
#MD 1 hlist_0.74980_SAM_Nb_0.fits
#MD 25 hlist_0.75440_SAM_Nb_10.fits

def create_plots(env='MD04'):
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "catalogs", "*.Ms.fits")))
	print fileList
	for fileN in fileList:
		print fileN
		hd = fits.open(fileN)[1].data	
		p.figure(1, (6,6))
		p.plot(hd['Mgal_mvir_Mo13'], hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'], 'b,', rasterized=True)
		p.xlabel(r'$\log_{10} M_\odot$')
		p.ylabel(r'$\log_{10}\lambda_{SAR} + \log_{10}M_\odot$')
		p.grid()
		p.savefig(os.path.join(os.environ[env], "results", os.path.basename(fileN)[:-5]+'.pdf'))
        p.clf()
        
create_plots(env='MD04')
create_plots(env='MD10')
create_plots(env='MD25')

os.system("cp $MD04/results/*.pdf ~/wwwDir/eRoMok/plots/MD_0.4Gpc")
os.system("cp $MD10/results/*.pdf ~/wwwDir/eRoMok/plots/MD_1.0Gpc")
os.system("cp $MD25/results/*.pdf ~/wwwDir/eRoMok/plots/MD_2.5Gpc")

