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

inFiles = n.array(["dmdens_cic_104_DFhist.dat", "dmdens_cic_101_DFhist.dat", "dmdens_cic_097_DFhist.dat", "dmdens_cic_087_DFhist.dat"])
# ZS = 0.7 0.8 1.0 1.48

def getNN0_sim(inSim,NR=10):
	f=open(join(mockDir, inSim))
	bins, HDF0 = cPickle.load(f)
	f.close()
	#bins = n.hstack((0,n.logspace(-3, 4, 1000)))
	xb = (bins[1:]+bins[:-1])/2.
	dx = bins[1:] - bins[:-1]
	X, Y = n.meshgrid(xb,xb)
	N0 = HDF0 /dx / (1000.-2*1000./2048)**3.
	HDF0R = n.array([HDF0[ii::NR] for ii in range(NR)]).sum(axis=0)
	binsR = bins[::NR]
	N0R = HDF0R / ((binsR[1:] - binsR[:-1]) * 250.**3.)
	return N0, bins, N0R, binsR


N0z07s, binsz07s, N0z07, binsz07 = getNN0_sim(inFiles[0])
xb = (binsz07[1:]+binsz07[:-1])/2.

f=open(join(mockDir,'Planck-ng512-L250.0.HDF0.pkl'),'r')
bins, HDF0, N0 = cPickle.load(f)
f.close()

NR = 10
HDF0R = n.array([HDF0[ii::NR] for ii in range(NR)]).sum(axis=0)
binsR = bins[::NR]
N0R = HDF0R / ((binsR[1:] - binsR[:-1]) * 250.**3.)
N0R_sig = n.array([N0[ii::NR] for ii in range(NR)]).std(axis=0)

muscleDelta = interp1d( N0R, (binsR[:-1]+binsR[1:])/2.)
mdplDelta = interp1d( N0z07, xb)
ok=(N0R>0)&(N0R<=100)#n.max(N0z07))
trueDelta = mdplDelta(N0R[ok])
index=n.argsort(N0R[ok])
deltaMuscle = (binsR[:-1]+binsR[1:])/2.
n.savetxt(join(mockDir,"delta-conversion-muscle-mdpl.txt"),n.transpose([deltaMuscle[ok][index],trueDelta[index]]),header="deltaMuscle deltaMDPL")

p.figure(0)
p.title('QSO')
p.plot(xb, N0z07,'kx', rasterized=True, label='z=0.7 all')
p.plot(xb, N0,'bx', rasterized=True, label='z=0.7 muscle')
p.plot((binsR[:-1]+binsR[1:])/2., N0R,'rx', rasterized=True, label='z=0.7 muscle resampled')
p.plot(trueDelta[index], N0R[ok][index], 'm--', lw=2, rasterized=True, label='z=0.7 muscle corr')
p.xlabel(r'$\delta_0$')
p.ylabel(r'N')
p.xscale('log')
p.yscale('log')
p.ylim((1e-10, 1e2))
p.xlim((0.1, 1e4))
gl = p.legend(loc=3)
gl.set_frame_on(False)
p.grid()
p.savefig(join(mockDir,"plots","muscle-delta-HDF0.png"))
p.show()