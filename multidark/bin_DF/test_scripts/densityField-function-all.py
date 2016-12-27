#  cd pySU/pyMultidark/trunk/bin/fortranfile-0.2.1/
import sys
import numpy as n
import os
from os.path import join
from astropy.io import fits
import time
import cPickle
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.stats import scoreatpercentile as sc
from scipy.stats import norm
from scipy.stats import lognorm
import matplotlib.pyplot as p
from matplotlib.ticker import NullFormatter
nullfmt = NullFormatter()         # no labels
from scipy.optimize import curve_fit
import astropy.cosmology as co
import astropy.units as uu

Mp = 1.51 * 10**9. * uu.solMass
CELL =1000./2048.*uu.megaparsec

DFdir = join("/data2", "users", "gustavo", "BigMD", "1Gpc_3840_Planck1_New", "DENSFIELDS")
# mockDir = "/data1/DATA/eBOSS/Multidark-box-mocks/parts/"
mockDir = join("..","MD_1Gpc","density_field")

def getNN0_sim(inSim):
	f=open(join(mockDir, inSim))
	bins, HDF0 = cPickle.load(f)
	f.close()
	#bins = n.hstack((0,n.logspace(-3, 4, 1000)))
	xb = (bins[1:]+bins[:-1])/2.
	dx = bins[1:] - bins[:-1]
	X, Y = n.meshgrid(xb,xb)
	N0 = HDF0 /dx #/ (1000.-2*1000./2048)**3.
	return N0, bins

def getNN0_sim_logBin(inSim):
	f=open(join(mockDir, inSim))
	bins, HDF0 = cPickle.load(f)
	f.close()
	#bins = n.hstack((0,n.logspace(-3, 4, 1000)))
	xb = (bins[1:]+bins[:-1])/2.
	dlogx = (bins[1:] - bins[:-1])/xb
	X, Y = n.meshgrid(xb,xb)
	N0 = HDF0 /dlogx #/ (1000.-2*1000./2048)**3.
	return N0, bins



cd = lambda z:  co.Planck13.Om(z) * co.Planck13.critical_density(z).to(uu.solMass/uu.megaparsec**3) * CELL**3./ (Mp * (co.Planck13.H(z)/100.)**2.)

inFiles = n.array(["dmdens_cic_104_DFhist_linearBins.dat", "dmdens_cic_101_DFhist_linearBins.dat", "dmdens_cic_097_DFhist_linearBins.dat", "dmdens_cic_087_DFhist_linearBins.dat"])

N0z07, binsz07 = getNN0_sim(inFiles[0])
bins = binsz07
xb07 = (bins[1:]+bins[:-1])/2. #* 100

N0z08, binsz08 = getNN0_sim(inFiles[1])
bins = binsz08
xb08 = (bins[1:]+bins[:-1])/2. #* 100

N0z15, binsz15 = getNN0_sim(inFiles[3])
bins = binsz15
xb15 = (bins[1:]+bins[:-1])/2.#* 100

dx = bins[1:] - bins[:-1]

dfn = lambda DF, DFs, a, A : n.log10( A * (10**DF/10**DFs)**(-a) *n.e**(-10**DF/10**DFs) )
ok = (N0z07[2:]*1e-9>0)
xfit = n.log10(xb07[2:][ok])
yfit = n.log10( N0z07[2:][ok]*1e-9 )
out07, cov = curve_fit(dfn, xfit, yfit, p0=(3, 2, 1e-6 ), maxfev=80000000)

ok = (N0z15[2:]*1e-9>0)
xfit = n.log10(xb15[2:][ok])
yfit = n.log10( N0z15[2:][ok]*1e-9 )
out15, cov = curve_fit(dfn, xfit, yfit, p0=(3, 2, 1e-6 ), maxfev=80000000)


p.figure(0)
p.plot(xb07[2:], N0z07[2:]*1e-9,'k', lw=2, rasterized=True, label='z=0.7 all')
p.plot(xb07[1], N0z07[1]*1e-9,'ko', rasterized=True, label='z=0.7 empty')
p.plot(xb07[2:], 10**dfn(n.log10(xb07[2:]), out07[0], out07[1], out07[2]), 'm--', label='fit z=0.7')
p.plot(xb08[2:], N0z08[2:]*1e-9,'b', lw=2, rasterized=True, label='z=0.8 all')
p.plot(xb08[1], N0z08[1]*1e-9,'bo', rasterized=True, label='z=0.8 empty')
p.plot(xb15[2:], N0z15[2:]*1e-9,'r',lw=2, rasterized=True, label='z=1.5 all')
p.plot(xb15[1], N0z15[1]*1e-9,'ro', rasterized=True, label='z=1.5 empty')
p.plot(xb15[2:], 10**dfn(n.log10(xb15[2:]), out15[0], out15[1], out15[2]), 'm--', label='fit z=1.5')
p.axvline(1/ cd(0.701838).value,label='1 particle per cell',c='r')
p.xlabel(r'density field value')
p.ylabel(r'N[in bin]/bin width/volume')
p.xscale('log')
p.yscale('log')
p.ylim((1e-10, 1e2))
#p.xlim((0.01, 1e4))
gl = p.legend(loc=0)
gl.set_frame_on(False)
p.title('cell of 0.48 Mpc/h')
p.grid()
p.savefig(join(mockDir,"plots","DF-HDF0-histogram_linear.png"))
p.clf()

#lgn = lambda x, mu, sig, A : A* n.e**(-(n.log(x)-mu)**2./(2*sig**2.)**0.5)/ (x * (2*sig**2.)**0.5)
#out, cov = curve_fit(lgn, xb15[(xb15>3)]* cd(1.480160).value, N0z15[(xb15>3)], p0=(7,1,1e9 ))

p.figure(0)
p.plot(xb07[2:] * cd(0.701838).value, N0z07[2:],'k', lw=2, rasterized=True, label='z=0.7 all')
p.plot(xb07[1] * cd(0.701838).value, N0z07[1],'ko', rasterized=True, label='z=0.7 empty')
p.plot(xb08[2:] * cd(0.818843).value , N0z08[2:],'b', lw=2, rasterized=True, label='z=0.8 all')
p.plot(xb08[1] * cd(0.818843).value , N0z08[1],'bo', rasterized=True, label='z=0.8 empty')
p.plot(xb15[2:] * cd(1.480160).value , N0z15[2:],'r',lw=2, rasterized=True, label='z=1.5 all')
p.plot(xb15[1] * cd(1.480160).value , N0z15[1],'ro', rasterized=True, label='z=1.5 empty')
#p.plot(xb15, lgn(xb15, out[0], out[1], out[2]), 'm--', label='lognormal')
#p.hist(pd,bins=bins, histtype='step')
p.axvline(1,label='1 particle per cell',c='r')
p.xlabel(r'N particles / cell')
p.ylabel(r'N[in bin]/bin width')
p.xscale('log')
p.yscale('log')
#p.ylim((0.9, 1e13))
#p.xlim((0.1, 1e4))
gl = p.legend(loc=0)
gl.set_frame_on(False)
p.title('cell of 0.48 Mpc/h')
p.grid()
p.savefig(join(mockDir,"plots","DF-Nparticle-histogram_linear.png"))
p.clf()

sys.exit()

inFiles = n.array(["dmdens_cic_104_DFhist.dat", "dmdens_cic_101_DFhist.dat", "dmdens_cic_097_DFhist.dat", "dmdens_cic_087_DFhist.dat"])


N0z07, binsz07 = getNN0_sim_logBin(inFiles[0])
bins = binsz07
xb07 = (bins[1:]+bins[:-1])/2. #* 100

N0z08, binsz08 = getNN0_sim_logBin(inFiles[1])
bins = binsz08
xb08 = (bins[1:]+bins[:-1])/2. #* 100

N0z15, binsz15 = getNN0_sim_logBin(inFiles[3])
bins = binsz15
xb15 = (bins[1:]+bins[:-1])/2.#* 100

dx = bins[1:] - bins[:-1]

p.figure(0)
p.plot(xb07[2:], N0z07[2:]*1e-9,'k', lw=2, rasterized=True, label='z=0.7 all')
p.plot(xb07[1], N0z07[1]*1e-9,'ko', rasterized=True, label='z=0.7 empty')
p.plot(xb08[2:], N0z08[2:]*1e-9,'b', lw=2, rasterized=True, label='z=0.8 all')
p.plot(xb08[1], N0z08[1]*1e-9,'bo', rasterized=True, label='z=0.8 empty')
p.plot(xb15[2:], N0z15[2:]*1e-9,'r',lw=2, rasterized=True, label='z=1.5 all')
p.plot(xb15[1], N0z15[1]*1e-9,'ro', rasterized=True, label='z=1.5 empty')
p.axvline(1/ cd(0.701838).value,label='1 particle per cell',c='r')
p.xlabel(r'density field value')
p.ylabel(r'dN/dlogbin/volume')
p.xscale('log')
p.yscale('log')
#p.ylim((0.9, 1e13))
#p.xlim((0.01, 1e4))
gl = p.legend(loc=0)
gl.set_frame_on(False)
p.title('cell of 0.48 Mpc/h')
p.grid()
p.savefig(join(mockDir,"plots","DF-HDF0-histogram_log.png"))
p.clf()

#lgn = lambda x, mu, sig, A : A* n.e**(-(n.log(x)-mu)**2./(2*sig**2.)**0.5)/ (x * (2*sig**2.)**0.5)
#out, cov = curve_fit(lgn, xb15[(xb15>3)]* cd(1.480160).value, N0z15[(xb15>3)], p0=(7,1,1e9 ))

p.figure(0)
p.plot(xb07[2:] * cd(0.701838).value, N0z07[2:],'k', lw=2, rasterized=True, label='z=0.7 all')
p.plot(xb07[1] * cd(0.701838).value, N0z07[1],'ko', rasterized=True, label='z=0.7 empty')
p.plot(xb08[2:] * cd(0.818843).value , N0z08[2:],'b', lw=2, rasterized=True, label='z=0.8 all')
p.plot(xb08[1] * cd(0.818843).value , N0z08[1],'bo', rasterized=True, label='z=0.8 empty')
p.plot(xb15[2:] * cd(1.480160).value , N0z15[2:],'r',lw=2, rasterized=True, label='z=1.5 all')
p.plot(xb15[1] * cd(1.480160).value , N0z15[1],'ro', rasterized=True, label='z=1.5 empty')
#p.plot(xb15, lgn(xb15, out[0], out[1], out[2]), 'm--', label='lognormal')
#p.hist(pd,bins=bins, histtype='step')
p.axvline(1,label='1 particle per cell',c='r')
p.xlabel(r'N particles / cell')
p.ylabel(r'dN/dlogbin')
p.xscale('log')
p.yscale('log')
#p.ylim((0.9, 1e13))
#p.xlim((0.1, 1e4))
gl = p.legend(loc=0)
gl.set_frame_on(False)
p.title('cell of 0.48 Mpc/h')
p.grid()
p.savefig(join(mockDir,"plots","DF-Nparticle-histogram_log.png"))
p.clf()




