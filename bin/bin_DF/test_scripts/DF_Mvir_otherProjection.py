
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
#DFdir = join("..", "MDPL")
DFfile = join("dmdens_cic_087.dat")
Halofile = join("hlist_0.40320_PM.DF.fits")
Lbox = 1000.
grid = 2048

qty = 'mvir'
dx = Lbox/grid


NNs =[8, 16, 32, 64 ]
massRanges=[[11,12],[12,13],[13,14]]
for NN in NNs:
	for massRange in massRanges:
		print NN, massRange
		massMin = massRange[0]
		massMax = massRange[1]
		dxN = NN * dx

		allSlices = n.array(glob.glob( "/users/jcomparat/skies_universes_production/MultiDark/onePtFct/slice"+str(NN)+"/occupationDATA_raw/*"+str(massMin)+".M."+str(massMax)+".pkl"))

		for jjj, path_to_slice in enumerate (allSlices):
			path_to_plot = path_to_slice.replace("occupationDATA_raw", "meanDelta_plot")[:-4]+".png"
			path_to_plot_d = path_to_slice.replace("occupationDATA_raw", "distrDelta_plot")[:-4]+".png"
			f=open(path_to_slice, 'r')
			DF_rs,NH = cPickle.load(f)
			f.close()
			#figure halos
			binsNH = n.arange(-0.5,n.max(NH)+1,1)
			xNH = (binsNH[1:]+binsNH[:-1])/2.
			deltaMass = n.zeros(len(binsNH[:-1]))
			deltaMassSTD = n.zeros(len(binsNH[:-1]))
			#compute occupation
			
			for numBin, bin in enumerate(binsNH[:-1]):
				ids = n.where( (NH>=binsNH[numBin])&(NH<binsNH[numBin+1]) )
				print bin, len(DF_rs[ids])
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
			p.savefig(path_to_plot)
			p.clf()
			
			p.figure(1)
			for numBin, bin in enumerate(binsNH[:-1]):
				ids = n.where( (NH>=binsNH[numBin])&(NH<binsNH[numBin+1]) )
				#print numBin, len(DF_rs[ids]), DF_rs[ids]
				if len(DF_rs[ids])>10 :
					print
					p.hist(DF_rs[ids], normed=True, label=str(xNH[numBin]), histtype='step', bins=20)

			
			p.legend(loc=0, fontsize=10)
			p.xlabel('log(1+delta)')
			p.ylabel('normed histogram')
			p.title("cell="+str(dxN)+"Mpc/h")
			p.grid()
			p.savefig(path_to_plot_d)
			p.clf()

