from os.path import join
import os
import glob
import time

import cPickle
import fileinput
import astropy.io.fits as fits

import numpy as n
from scipy.interpolate import interp1d
import scipy.spatial.ckdtree as t

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

plotDir = os.path.join(os.environ['HOME'], 'wwwDir', "eRoMok", "snapshot_sanity_checks")

sf = fits.open(os.path.join(os.environ["MD10"], "output_MD_1.0Gpc.fits"))[1].data

for name in sf['snap_name']:
	path_2_file = os.path.join(os.environ["MD10"], 'work_agn', "out_" + name + "_SAM_Nb_0.fits")
	print path_2_file
	if os.path.isfile(path_2_file):
		hdus = fits.open(path_2_file)
		for qty_name in hdus[1].columns.names[:-3] :
			f=open(os.path.join(os.environ["MD10"], 'work_agn', 'sanity_check', os.path.basename(path_2_file)[:-5] + "_" + qty_name + ".pkl"),'r')
			out = cPickle.load(f)
			f.close()
			
			p.figure(1, (6,6))
			p.plot((out[1][:-1]+out[1][1:])/2. , n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ])[::-1] )
			p.xlabel(qty_name)
			p.ylabel('counts')
			p.yscale('log')
			p.grid()
			p.savefig(os.path.join(plotDir, os.path.basename(path_2_file)[:-5] + "_" + qty_name + ".jpg"))
			p.clf()
			
