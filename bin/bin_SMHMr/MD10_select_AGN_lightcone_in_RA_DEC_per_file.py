# overall python packages
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys
t0=time.time()


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
		agn = (data[3] > 30 ) & (data[3] < 40 ) & (data[4] > 30 ) & (data[4] < 40 )
		print len(agn.nonzero()[0])
		#data_all = n.vstack((data_all, data.T[agn]))
		#print data_all.shape
		n.savetxt(os.path.join(os.environ[env], "light-cone", os.path.basename(fileName)[:-4]+".small-agn-catalog.ascii"), data.T[agn])	

get_AGN_catalog(env='MD10')
print time.time()-t0, "seconds"


