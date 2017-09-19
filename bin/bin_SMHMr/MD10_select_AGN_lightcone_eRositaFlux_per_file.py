# overall python packages
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys
t0=time.time()
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)


def get_AGN_catalog(env='MD10'):
	# gets the file list to add the Xray luminosity
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "light-cone", "MDPL2_ROCKSTAR_FluxProj_*_000_AGN.dat" )))
	fileList.sort()
	print fileList
	#print fileList[0]
	#data = n.loadtxt(fileList[0],unpack=True)
	#print data, data.shape
	#agn = (data[3] > 30 ) & (data[3] < 40 ) & (data[4] > 30 ) & (data[4] < 40 )
	#data_all = data.T[agn]
	#print data_all.shape
	for fileName in fileList:
		print fileName
		data = n.loadtxt(fileName,unpack=True)
		print data.shape
		# compute luminosity
		dL_cm = (cosmoMD.luminosity_distance(data[2]).to(u.cm)).value
		flux = 10**(data[9]-0.3) / (4.* n.pi * dL_cm * dL_cm)
		print dL_cm, flux
		agn = (flux > 1e-15 ) #& (data[2] < 2.44)
		print len(agn.nonzero()[0])
		#data_all = n.vstack((data_all, data.T[agn]))
		#print data_all.shape
		n.savetxt(os.path.join(os.environ[env], "light-cone", os.path.basename(fileName)[:-4]+".erosita-agn-full-depth-catalog.ascii"), data.T[agn])	

get_AGN_catalog(env='MD10')
print time.time()-t0, "seconds"


