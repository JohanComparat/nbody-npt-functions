import StellarMass
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

plotDir = os.path.join(os.environ['HOME'], 'wwwDir', "eRoMok", "snapshot_sanity_checks")

def measureSMF(sf, env='MD10', volume=1000.**3., out_dir="../"):
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", "out_"+el["snap_name"]+"_SAM_Nb_*_Ms.fits")))
	fileList.sort()
	snapList = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", "out_"+el["snap_name"]+"_SAM_Nb_?.fits")))
	snapList.sort()
	print fileList, snapList
	for fileMs, fileN in zip(fileList, snapList):
		print fileN, fileMs
		hM = fits.open(fileMs)
		hN = fits.open(fileN)
		mass = hM[1].data['stellar_mass_Mo13_mvir']
		mvir = hN[1].data['mvir']
		smhmr = mass - mvir
		selection = (hM[1].data['stellar_mass_reliable'])
		
		p.figure(1, (6,6))
		p.plot(mvir[selection], smhmr[selection], 'k,', rasterized = True )
		p.plot(mvir[selection], smhmr[selection]+0.15/2., 'r,', rasterized = True )
		p.plot(mvir[selection], smhmr[selection]-0.15/2., 'r,', rasterized = True )
		p.xlabel('mvir')
		p.ylabel('stellar mass - mvir')
		#p.yscale('log')
		p.title(str(el['redshift']))
		p.grid()
		p.savefig(os.path.join(plotDir, os.path.basename(fileMs)[:-5] + "_smhmr_moster2013.jpg"))
		p.clf()


	
# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

out_dir = os.path.join(os.path.join(os.environ['MD10'], "results", "stellar_mass_function", "data"))

for el in summ:
	print el
	measureSMF(el, env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S5_BCG", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S5_GAL", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S6_AGN", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S8_BG1", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S8_BG2", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S8_ELG", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S8_QSO", env='MD10', volume=1000.**3., out_dir = out_dir)
	

