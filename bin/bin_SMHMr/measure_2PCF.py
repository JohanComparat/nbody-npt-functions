from MultiDark import *

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

box = MultiDarkSimulation(Lbox=1000.0 * u.Mpc, boxDir = "MD_1Gpc")

tracer_names = n.array(["4MOST_S5_BCG", "4MOST_S5_GAL", "4MOST_S6_AGN", "4MOST_S8_LRG", "4MOST_S8_ELG", "4MOST_S8_QSO"])

summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

out_dir = os.path.join(os.path.join(os.environ['MD10'],"results","mvir_mass_function", "data"))

def get_xyz(snap_name, tracer_name, env='MD10'):
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", "out_"+snap_name+"_SAM_Nb_?.fits")))
	fileList.sort()
	fileList_T = n.array(glob.glob(os.path.join(os.environ[env], "work_agn", "out_"+snap_name+"_SAM_Nb_*_"+tracer_name+".fits")))
	fileList_T.sort()
	tracer_name
	print fileList, fileList_T
	if len(fileList_T)==len(fileList):
		x, y, z = [], [], []
		for ii, fileN in enumerate(fileList):
			print fileN
			hh = fits.open(fileN)
			lines = fits.open(fileList_T[ii])[1].data['line_number']
			x.append( hh[1].data['x'][lines] )
			y.append( hh[1].data['y'][lines] )
			z.append( hh[1].data['z'][lines] )
		return n.hstack((x)), n.hstack((y)), n.hstack((z))
	else:
		return [-1], [-1], [-1]

for el in summ:
	print el
	out_file = lambda tracer_name : os.path.join(out_dir, "out_"+el['snap_name']+"_"+tracer_name+"_2PCF.pkl")
	for tracer_name in tracer_names :
		print tracer_name
		fileList_T = n.array(glob.glob(os.path.join(os.environ['MD10'], "work_agn", "out_"+el['snap_name']+"_SAM_Nb_*_"+tracer_name+".fits")))
		print len(fileList_T)
		if os.path.isfile(out_file(tracer_name))==False and len(fileList_T)>0:
			xR, yR, zR = get_xyz(el['snap_name'], tracer_name)
			print len(xR)
			if len(xR)>5000. and len(xR)<1000000. :
				box.compute2PCF_single(xR, yR, zR, outfile=out_file(tracer_name), rmax=30, dr = 0.5)
			elif len(xR)>1000000.:
				downSamp = (n.random.random(len(xR))<1000000. / float(len(xR)) )
				xR = xR[downSamp]
				yR = yR[downSamp]
				zR = zR[downSamp]
				box.compute2PCF_single(xR, yR, zR, outfile=out_file(tracer_name), rmax=30, dr = 0.5)
			else :
				pass



