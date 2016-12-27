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
import matplotlib.pyplot as p
# loads the density field
# listof parameters :
#DFdir = join("..", "MDPL")
DFfile = join("..", "MDPL", "dmdens_cic_087.dat")
Halofile = join("..", "MDPL", "hlist_0.40320_PM.DF.fits")

qty = 'mvir'
NN=32
Lbox = 1000.
grid = 2048
n_slices = grid/NN 
binsQTY = n.arange(11,14.7,0.1)
binsDelta = 10**n.arange(-0.2,0.2,0.002)
logbinsDelta = n.arange(-0.2,0.2,0.002)
dx = Lbox/grid
dxN = NN * dx

path_to_outputDF = Halofile[:-4] + ".degrN."+str(NN)+".pkl"
path_to_outputHIST = Halofile[:-7] + qty + "HIST.degrN."+str(NN)+".pkl"
path_to_outputMatrix = Halofile[:-4] + qty + ".degrN."+str(NN)+".matrix"

deltaMass = n.loadtxt(path_to_outputMatrix)
deltaMassSTD = n.loadtxt(path_to_outputMatrix+"std")

xs = (logbinsDelta[1:] + logbinsDelta[:-1])/2.
ys = (binsQTY[1:] + binsQTY[:-1])/2.
xlogMvir, xlogDelta = n.meshgrid(ys, xs)

#p.contour(xlogMvir, xlogDelta,n.log10(1+deltaMass))
#p.show()

test = n.sum(deltaMass, axis=1)
test_std = n.sum(deltaMassSTD, axis=1)

p.plot(ys, deltaMass[-10])
p.plot(ys, deltaMass[10])

p.plot(10**xs, test)
p.plot(10**xs, test-test_std, 'r--')
p.plot(10**xs, test+test_std, 'r--')
p.xlabel('1+delta DM')
p.ylabel('average Nha (log Mvir>11)')
p.show()


p.figure(0,(6,6))
p.scatter( n.ravel(xlogMvir), 1+n.ravel(deltaMass), c=n.ravel(xlogDelta), s =5, edgecolors='face' )
cb = p.colorbar()
cb.set_label('log(1+delta DM)')
p.xlabel('log Mvir')
p.ylabel('1+N halo')
p.yscale('log')
p.ylim(0.9,1e5)
p.title("cell="+str(n.round(dxN,2))+"Mpc/h")
p.show()

p.figure(0,(6,6))
p.scatter( n.ravel(xlogDelta), 1+n.ravel(deltaMass), c=n.ravel(xlogMvir), s =5, edgecolors='face' )
cb = p.colorbar()
cb.set_label('log Mvir')
p.xlabel('log(1+delta DM)')
p.ylabel('1+N halo')
p.yscale('log')
p.ylim(0.9,1e5)
p.title("cell="+str(n.round(dxN,2))+"Mpc/h")
p.show()



