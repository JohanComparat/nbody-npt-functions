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

sf = fits.open(os.path.join(os.environ["MD10"], "output_MD_1.0Gpc.fits"))[1].data

for name in sf['snap_name']:
	path_2_file = os.path.join(os.environ["MD10"], 'work_agn', "out_" + name + "_SAM_Nb_0_Ms.fits")
	print path_2_file
	if os.path.isfile(path_2_file):
		hdus = fits.open(path_2_file)
		ok = (hdus[1].data['stellar_mass_reliable'])

		for qty_name in hdus[1].columns.names[:-1] :
			t0 = time.time()
			array_qty = hdus[1].data[qty_name][ok]
			print qty_name, n.min(array_qty), n.max(array_qty)
			out = n.histogram(array_qty, bins = 1000)
			f=open(os.path.join(os.environ["MD10"], 'work_agn', 'sanity_check', os.path.basename(path_2_file)[:-5] + "_" + qty_name + ".pkl"),'w')
			cPickle.dump(out, f)
			f.close()
			print "time elapsed", time.time() - t0
