#  cd pySU/pyMultidark/trunk/bin/fortranfile-0.2.1/

import numpy as n
import os
from os.path import join
from astropy.io import fits
import time
import fortranfile
import glob
import cPickle
import sys
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
# loads the density field
# listof parameters :
DF_data_dir = join(os.environ['MD_DF_DIR'], "data")
DF_results_dir = join(os.environ['MD_DF_DIR'], "results")

DFfile = join(DF_data_dir, "dmdens_cic_087.dat")
Halofile = join(DF_data_dir, "hlist_0.40320_PM.DF.fits")

Lbox = 1000.
grid = 2048

qty = 'mvir'
dx = Lbox/grid


#binsQTY = n.arange(11,14.7,0.1)
def computeSlice(NN=4, massRange=[11,12], sliceNum=0, 	logbinsDelta = n.arange(-1., 2.6, 0.1), cos='cen_sat'):
	massMin = massRange[0]
	massMax = massRange[1]
	dxN = NN * dx
	n_slices_total = grid/NN 
	n_slices = 1 # grid/NN 

	path_to_outputSlice = join(DF_results_dir, "slice"+str(NN), "occupationDATA_raw", Halofile.split('/')[-1][:-8] + ".slice"+str(sliceNum)+"."+str(massMin)+".M."+str(massMax)+".pkl")
	path_to_outputMatrix = join(DF_results_dir, "slice"+str(NN), "occupationDATA", Halofile.split('/')[-1][:-8] + ".slice"+str(sliceNum)+"."+str(massMin)+".M."+str(massMax)+".matrix")
	
	#open halo file
	md = fits.open(Halofile)[1].data
	# function 1
	def getSliceNhalo(NN=32, sliceNum=0):
		dx = Lbox/grid
		dxN = NN * dx
		ix  = ( ( md['x'] / dxN ) // 1 ).astype( 'int' )
		iy  = ( ( md['y'] / dxN ) // 1 ).astype( 'int' )
		iz = ( ( md['z'] / dxN ) // 1 ).astype( 'int' )
		if cos=='cen':
			sel = (iz==sliceNum)&(md['pid']<0)&(md[qty]>massMin)&(md[qty]<massMax)
		if cos=='cen_sat':
			sel = (iz==sliceNum)&(md[qty]>massMin)&(md[qty]<massMax)
		H, xedges, yedges = n.histogram2d(ix[sel], iy[sel], bins=n.arange(0,grid/NN+1,1))
		return H.T
	# function 2
	def getSliceDF(NN=32, sliceNum=0):
		if sliceNum==0:
			f = fortranfile.FortranFile(DFfile)
			gridx, gridy, gridz = f.readInts()
			#DF_sl1 = n.zeros((NN, gridx, gridy))
			DF_sl1 = n.array([f.readReals().reshape(2048,2048) for ii in range(NN)])
			return n.mean(DF_sl1, axis=0)
		if sliceNum>0:
			f = fortranfile.FortranFile(DFfile)
			gridx, gridy, gridz = f.readInts()
			for sli in n.arange(sliceNum):
				for ii in range(NN):
					f.readReals()
			DF_sl1 = n.array([f.readReals().reshape(2048,2048) for ii in range(NN)])
			return n.mean(DF_sl1, axis=0)

	# get 2d mesh
	DF = getSliceDF(NN, sliceNum)
	#create x, y mesh
	jjMesh,kkMesh = n.meshgrid(n.arange(grid/NN), n.arange(grid/NN))
	js = n.ravel(jjMesh)
	ks = n.ravel(kkMesh)
	# rebins mesh
	DF_rs = n.array([ n.mean(DF[jj*NN:jj*NN+NN:1, kk*NN:kk*NN+NN:1]) for jj,kk in zip(js, ks)]).reshape(grid/NN, grid/NN).T
	#gets halo slice
	NH = getSliceNhalo(NN, sliceNum)
	#writes results
	f=open(path_to_outputSlice, 'w')
	cPickle.dump([DF_rs,NH], f)
	f.close()
	f=open(path_to_outputSlice, 'r')
	DF_rs,NH = cPickle.load(f)
	f.close()
	#figure halos
	p.figure(1)
	p.imshow(n.log10(0.0001+NH), vmin=-0.5, vmax=1.5)
	p.colorbar()
	p.title("log10(Nhalo)")
	p.savefig(join(DF_results_dir, "slice"+str(NN), "Nhalo_plot", Halofile.split('/')[-1][:-8] + ".Nhalo.slice"+str(sliceNum)+"."+str(massMin)+".M."+str(massMax)+".png"))
	p.clf()
	#figure delta
	p.figure(1)
	p.imshow(n.log10(0.0001+DF_rs), vmin=-0.5, vmax=1.5)
	p.colorbar()
	p.title("log10(1+delta DM)")
	p.savefig(join(DF_results_dir, "slice"+str(NN), "DF_plot", Halofile.split('/')[-1][:-8] + ".DF.slice"+str(sliceNum)+"."+str(massMin)+".M."+str(massMax)+".png"))
	p.clf()
	
	#compute occupation
	xDelta = (logbinsDelta[1:]+logbinsDelta[:-1])/2.
	deltaMass = n.zeros(len(logbinsDelta[:-1]))
	deltaMassSTD = n.zeros(len(logbinsDelta[:-1]))
	for numBin, bin in enumerate(logbinsDelta[:-1]):
		ids = n.where( (DF_rs>=10**logbinsDelta[numBin])&(DF_rs<10**logbinsDelta[numBin+1]) )
		#print n.mean(NH[ids], axis=0), n.sum(NH[ids], axis=0)/float(len(NH[ids]))
		deltaMass[numBin] = n.mean(NH[ids], axis=0)
		deltaMassSTD[numBin] = n.std(NH[ids], axis=0)
	# figure occupation
	p.figure(1)
	p.plot(xDelta, deltaMass, 'b')
	p.plot(xDelta, deltaMass-deltaMassSTD, 'b--')
	p.plot(xDelta, deltaMass+deltaMassSTD, 'b--')
	p.xlabel('log10(1+delta)')
	p.ylabel('mean number of halos with '+str(massMin)+".M."+str(massMax))
	p.title("cell="+str(dxN)+"Mpc/h")
	p.grid()
	p.savefig(join(DF_results_dir, "slice"+str(NN), "meanOccupation_plot", Halofile.split('/')[-1][:-8] + ".meanOccupation.slice"+str(sliceNum)+"."+str(massMin)+".M."+str(massMax)+".png"))
	p.clf()
	# save occupation statistics
	n.savetxt(path_to_outputMatrix, n.transpose([xDelta,deltaMass, deltaMassSTD]), fmt='%10.5f', header="log10(1+delta) meanOccupation STDdev")

	#figure halos
	binsNH = n.arange(-0.5,n.max(NH)+1,1)
	xNH = (binsNH[1:]+binsNH[:-1])/2.
	deltaMass = n.zeros(len(binsNH[:-1]))
	deltaMassSTD = n.zeros(len(binsNH[:-1]))
	#compute occupation
	
	for numBin, bin in enumerate(binsNH[:-1]):
		ids = n.where( (NH>=binsNH[numBin])&(NH<binsNH[numBin+1]) )
		#print bin, len(DF_rs[ids])
		deltaMass[numBin] = n.mean(DF_rs[ids], axis=0)
		deltaMassSTD[numBin] = n.std(DF_rs[ids], axis=0)

	# figure occupation
	p.figure(1)
	p.plot(xNH, deltaMass, 'b')
	p.plot(xNH, deltaMass-deltaMassSTD, 'b--')
	p.plot(xNH, deltaMass+deltaMassSTD, 'b--')
	p.xlabel('N halo per cell in '+str(massMin)+".M."+str(massMax))
	p.ylabel('mean value of DF')
	p.title("cell="+str(dxN)+"Mpc/h")
	p.grid()
	p.savefig(join(DF_results_dir, "slice"+str(NN), "meanDelta_plot", Halofile.split('/')[-1][:-8] + ".meanDelta.slice"+str(sliceNum)+"."+str(massMin)+".M."+str(massMax)+".png"))
	p.clf()

	p.figure(1)
	for numBin, bin in enumerate(binsNH[:-1]):
		ids = n.where( (NH>=binsNH[numBin])&(NH<binsNH[numBin+1]) )
		#print numBin, len(DF_rs[ids]), DF_rs[ids]
		if len(DF_rs[ids])>10 :
			#print
			p.hist(DF_rs[ids], normed=True, label=str(xNH[numBin]), histtype='step', bins=20)

	
	p.legend(loc=0, fontsize=10)
	p.xlabel('log(1+delta)')
	p.ylabel('normed histogram')
	p.title("cell="+str(dxN)+"Mpc/h")
	p.grid()
	p.savefig(join(DF_results_dir, "slice"+str(NN), "distrDelta_plot", Halofile.split('/')[-1][:-8] + ".distrDelta.slice"+str(sliceNum)+"."+str(massMin)+".M."+str(massMax)+".png"))
	p.clf()


	
	
		
NNs =[8 ]
sliceNums=n.arange(256)
massRanges=[[11,12],[12,13],[13,14]]
for NN in NNs:
	for sliceNum in sliceNums:
		for massRange in massRanges:
			print NN, massRange, sliceNum
			computeSlice(NN, massRange, sliceNum)
