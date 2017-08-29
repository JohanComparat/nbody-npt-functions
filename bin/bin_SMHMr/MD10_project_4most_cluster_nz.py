import numpy as n
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d

import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

mbins = n.arange(8,14.5,0.25)

summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

out_dir = os.path.join(os.path.join(os.environ['MD10'],"results","mvir_mass_function", "data"))
data = []
for el in summ[2:]:
	#print el
	fileList_G = n.array(glob.glob(os.path.join(os.environ['MD10'], "work_agn", "out_"+el['snap_name']+"_SAM_Nb_*_4MOST_S5_GAL.fits")))
	fileList_B = n.array(glob.glob(os.path.join(os.environ['MD10'], "work_agn", "out_"+el['snap_name']+"_SAM_Nb_*_4MOST_S5_BCG.fits")))
	NN=[]
	for f1, f2 in zip(fileList_G, fileList_B):
		NN.append( [len(fits.open(f1)[1].data), len(fits.open(f2)[1].data)] )
		
	out = n.array( NN ) 
	#print out
	total = n.sum(out, axis=0)
	#print total
	print el['redshift'], 1.*total[0]/total[1]
	data.append([el['redshift'], total[0], total[1]])

etc = fits.open(os.path.join(os.environ['HOME'], 'data-4most/comparat/sim/round9/d/run01/ETC_20min/', 'output_catalogue_S5_LR.fits'))[1].data


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

path_2_fig = os.path.join(os.environ['HOME'], 'wwwDir', '4FS_texp', 'cluster_texp.png')

p.figure(0,(5.5,5.5))
p.hist(etc['TEXP_D'], bins=n.arange(0,200,1), histtype='step', normed=True, cumulative=True, label='Dark time')
p.hist(etc['TEXP_B'][etc['TEXP_B']>0], bins=n.arange(0,200,1), histtype='step', normed=True, cumulative=True, label='Bright time')
p.axvline(20, label='20 min', color='k')
p.axvline(120, label='2 h', color='r')
p.xscale('log')
p.ylabel('Normed cumulative histogram')
p.xlabel('exposure time [minutes]')
p.xlim((0.9,200))
p.ylim((-0.01,1.01))
p.grid()
p.title('S5 clusters')
p.legend(frameon=False, loc=0)
p.savefig(path_2_fig)
p.clf()

path_2_fig = os.path.join(os.environ['HOME'], 'wwwDir', '4FS_texp', 'cluster_mags.png')

p.figure(0,(5.5,5.5))
p.hist(etc['R_MAG'], bins=n.arange(17.2,22.5,0.1), histtype='step', normed=True, cumulative=True)
p.ylabel('Normed cumulative histogram')
p.xlabel('SDSS R magnitude')
p.xlim((17.2,22.5))
p.ylim((-0.01,1.01))
p.grid()
p.title('S5 clusters')
p.legend(frameon=False, loc=0)
p.savefig(path_2_fig)
p.clf()

