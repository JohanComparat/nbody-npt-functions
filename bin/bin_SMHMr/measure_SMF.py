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

#duty_cycle = 0.01

bins = n.arange(6,13,0.1)
xb = (bins[1:] + bins[:-1]) / 2.

def measureSMF(snap_name, env='MD10', volume=1000.**3., out_dir="../"):
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", "out_"+snap_name+"_SAM_Nb_*_Ms.fits")))
	fileList.sort()
	print fileList
	Hall = n.zeros((len(fileList), len(bins)-1))
	for ii, fileN in enumerate(fileList):
		print fileN
		hh = fits.open(fileN)
		mass = hh[1].data['stellar_mass_Mo13_mvir']
		selection = (mass>0) # (hh[1].data['stellar_mass_reliable'])&(mass>0)
		Hall[ii], bb = n.histogram(hh[1].data['stellar_mass_Mo13_mvir'], bins=bins)
	
	counts = n.sum(Hall, axis=0)
	dN_dVdlogM = counts*0.6777**3./(bins[1:]-bins[:-1])/volume/n.log(10)
	data = n.transpose([bins[:-1], bins[1:], counts, dN_dVdlogM ])
	n.savetxt(os.path.join(out_dir, "out_"+snap_name+"_SMF.txt"), data, header = "logMs_low logMs_up counts dN_dVdlogM")


def measureSMF_tracer(snap_name, tracer_name, env='MD10', volume=1000.**3., out_dir="../"):
		out_file = os.path.join(out_dir, "out_"+snap_name+"_"+tracer_name+"_SMF.txt")
		#if os.path.isfile(out_file)==False:
		fileList = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", "out_"+snap_name+"_SAM_Nb_*_Ms.fits")))
		fileList.sort()
		fileList_T = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", "out_"+snap_name+"_SAM_Nb_*_"+tracer_name+".fits")))
		fileList_T.sort()
		tracer_name
		print fileList, fileList_T
		if len(fileList_T)==len(fileList):
			Hall = n.zeros((len(fileList), len(bins)-1))
			for ii, fileN in enumerate(fileList):
				print fileN
				hh = fits.open(fileN)
				lines = fits.open(fileList_T[ii])[1].data['line_number']
				mass = hh[1].data['stellar_mass_Mo13_mvir'][lines]
				Hall[ii], bb = n.histogram(mass, bins=bins)
			
			counts = n.sum(Hall, axis=0)
			dN_dVdlogM = counts*0.6777**3./(bins[1:]-bins[:-1])/volume/n.log(10)
			data = n.transpose([bins[:-1], bins[1:], counts, dN_dVdlogM ])
			n.savetxt(out_file, data, header = "logMs_low logMs_up counts dN_dVdlogM")


# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

out_dir = os.path.join(os.path.join(os.environ['MD10'], "results", "stellar_mass_function", "data"))

for el in summ:
	print el
	measureSMF(snap_name=el["snap_name"], env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S5_BCG", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S5_GAL", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S6_AGN", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S8_BG1", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S8_BG2", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S8_ELG", env='MD10', volume=1000.**3., out_dir = out_dir)
	#measureSMF_tracer(snap_name=el["snap_name"], tracer_name="4MOST_S8_QSO", env='MD10', volume=1000.**3., out_dir = out_dir)
	

