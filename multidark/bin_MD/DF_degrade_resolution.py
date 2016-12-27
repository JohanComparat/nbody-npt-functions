#  cd pySU/pyMultidark/trunk/bin/fortranfile-0.2.1/

import numpy as n
import os
from os.path import join
from astropy.io import fits
import time
import fortranfile
import cPickle
DFdir = join("/data2", "users", "gustavo", "BigMD", "1Gpc_3840_Planck1_New", "DENSFIELDS")
mockDir = "/data1/DATA/eBOSS/Multidark-box-mocks/v1.0/parts/"
inFiles = n.array(["dmdens_cic_104.dat", "dmdens_cic_101.dat", "dmdens_cic_097.dat", "dmdens_cic_087.dat"])

bins = n.hstack((0,n.logspace(-3, 4, 1000)))

for infi in inFiles:
	print infi
	DFfile = join(DFdir,infi)
	f = fortranfile.FortranFile(DFfile)
	gridx, gridy, gridz = f.readInts()
	Ntot = gridx/2
	res0 = n.empty((Ntot, len(bins)-1))
	NS = n.arange(Ntot)
	for kk in NS:
		print kk, time.time()
		DFa = f.readReals()
		DFb = f.readReals()
		DFaR = DFa.reshape((gridx, gridx))
		DFbR = DFb.reshape((gridx, gridx))
		DF = n.mean(n.array([DFaR,DFbR]), axis=0)
		DFdg = n.array([ n.array([ n.mean([DF[2*i][2*j:2*j+2], DF[2*i+1][2*j:2*j+2]]) for j in NS]) for i in NS])
		res0[kk] = n.histogram(n.hstack((DFdg)), bins=bins)[0]

	f.close()
	path_to_outputCat =  join(mockDir,infi[:-4] + "_DF_dg2_hist.dat")
	f=open(path_to_outputCat, 'w')
	cPickle.dump( [bins, n.sum(res0, axis=0)], f )
	f.close()

	
sys.exit()

inFiles = n.array(["vx_cic_104.dat", "vx_cic_101.dat", "dmdens_cic_097.dat", "dmdens_cic_087.dat"])

inFiles = n.array(["vx_cic_104.dat", "vx_cic_101.dat", "vx_cic_097.dat", "vx_cic_087.dat", "vy_cic_104.dat", "vy_cic_101.dat", "vy_cic_097.dat", "vy_cic_087.dat", "vz_cic_104.dat", "vz_cic_101.dat", "vz_cic_097.dat", "vz_cic_087.dat"])
bins = n.arange(-2000.,2000., 5.)


for infi in inFiles:
	print infi
	DFfile = join(DFdir,infi)
	f = fortranfile.FortranFile(DFfile)
	gridx, gridy, gridz = f.readInts()
	res0 = n.empty((gridx, len(bins)-1))
	res1 = n.empty((gridx, len(bins)-1))
	resH = n.empty((gridx, len(bins)-1, len(bins)-1))
	for kk in range(gridx):
		DF = f.readReals()
		i = n.arange(1, gridx-1, 1)
		j = n.arange(1, gridx-1, 1)
		DF0 = DF[n.hstack((n.outer(i,j)))]
		N1 = n.transpose([ n.hstack((n.outer(i-1,j-1))), n.hstack((n.outer(i,j-1))), n.hstack((n.outer(i-1,j))), n.hstack((n.outer(i+1,j+1))), n.hstack((n.outer(i+1,j))), n.hstack((n.outer(i,j+1))), n.hstack((n.outer(i+1,j+1))), n.hstack((n.outer(i-1,j+1))) ])
		#  N1 = n.transpose([ (i-1) + gridx * (j -1), (i) + gridx * (j -1), (i-1) + gridx * (j), (i+1) + gridx * (j +1), (i+1) + gridx * (j ), (i) + gridx * (j +1), (i+1) + gridx * (j -1), (i-1) + gridx * (j +1) ]) 
		DF1 = n.array([ n.mean(DF[el]) for el in N1 ]) 
		res0[kk] = n.histogram(DF0,bins=bins)[0]
		res1[kk] = n.histogram(DF1,bins=bins)[0]
		resH[kk] = n.histogram2d(DF0, DF1, bins)[0]

	f.close()
	path_to_outputCat =  join(mockDir,infi[:-4] + "_DF0DF1hist.dat")
	f=open(path_to_outputCat, 'w')
	cPickle.dump([bins,n.sum(res0,axis=0), n.sum(res1,axis=0), n.sum(resH,axis=0)],f)
	f.close()
