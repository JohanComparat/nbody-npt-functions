import cPickle
import numpy as n
import astropy.cosmology as co
import astropy.units as uu
aa =co.Planck13
import time
from astropy.io import fits
import os 
from os.path import join
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
mockDir = join("..","MD_1Gpc","density_field")

fi1 = join(mockDir, "dmdens_cic_104_DFhist.dat")
fi2 = join(mockDir, "parts", "dmdens_cic_104_DF_dg4_hist.dat")
f=open(join(mockDir,'Planck-ng512-L250.0.HDF0.pkl'),'r')
binsMU, HDF0MU, N0MU = cPickle.load(f)
xbMU = (binsMU[1:]+binsMU[:-1])/2.
f.close()

# ZS = 0.7 

def getNN0_sim(file,NR=10):
	f=open(file)
	bins, HDF0 = cPickle.load(f)
	f.close()
	#bins = n.hstack((0,n.logspace(-3, 4, 1000)))
	xb = (bins[1:]+bins[:-1])/2.
	dx = bins[1:] - bins[:-1]
	X, Y = n.meshgrid(xb,xb)
	N0 = HDF0 /dx / (1000.)**3.
	return N0, bins

N0z07s, binsz07s = getNN0_sim(fi1)
xb = (binsz07s[1:]+binsz07s[:-1])/2.

f=open(fi2,'r')
binsDG, histVal = cPickle.load(f)
f.close()
N0z07DG = histVal / ((binsDG[1:] - binsDG[:-1]) * (1000.)**3.)
xb2 = (binsDG[1:]+binsDG[:-1])/2.

p.figure(0)
p.title('MDPL meshes')
p.plot(xb, N0z07s,'kx', rasterized=True, label='z=0.7 2048')
p.plot(xb2, N0z07DG*4**3,'bx', rasterized=True, label='z=0.7 512')
#p.plot(xbMU, HDF0MU/((binsMU[1:]-binsMU[:-1])*512**3.),'rx', rasterized=True, label='z=0.7 Muscle 512')
p.xlabel(r'$\delta_0$')
p.ylabel(r'N')
p.xscale('log')
p.yscale('log')
p.ylim((1e-10, 1e2))
p.xlim((0.1, 1e4))
gl = p.legend(loc=3)
gl.set_frame_on(False)
p.grid()
p.savefig(join(mockDir,"plots","MDPL-delta-degraded.png"))
p.show()