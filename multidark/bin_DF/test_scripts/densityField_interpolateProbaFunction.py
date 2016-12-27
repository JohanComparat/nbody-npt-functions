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
import matplotlib.pyplot as p
from matplotlib.ticker import NullFormatter
nullfmt = NullFormatter()         # no labels
from scipy.optimize import curve_fit
DFdir = join("/data2", "users", "gustavo", "BigMD", "1Gpc_3840_Planck1_New", "DENSFIELDS")
# mockDir = "/data1/DATA/eBOSS/Multidark-box-mocks/parts/"
mockDir = join("..","MD_1Gpc","density_field")

#inFiles = n.array(["dmdens_cic_104_DFhist.dat",cd "dmdens_cic_101_DFhist.dat", "dmdens_cic_097_DFhist.dat", "dmdens_cic_087_DFhist.dat"])

# inFiles = n.array(["dmdens_cic_104_DF0DF1hist.dat", "dmdens_cic_101_DF0DF1hist.dat", "dmdens_cic_097_DF0DF1hist.dat", "dmdens_cic_087_DF0DF1hist.dat"])

inFiles = n.array(["dmdens_cic_104_DFhist.dat", "dmdens_cic_101_DFhist.dat", "dmdens_cic_097_DFhist.dat", "dmdens_cic_087_DFhist.dat"])
# ZS = 0.7 0.8 1.0 1.48

bins = n.logspace(-1, 4, 10000)
dx = bins[1:] - bins[:-1]
xb = (bins[1:]+bins[:-1])/2.
sel=(xb>0.1)&(xb<3000)
x = n.log10(xb[sel])

"z0.701838"
"z0.818843"
"z0.987281"
"z1.480160"

ps = n.loadtxt("fit-polynomial-ELG-z07.data")
proba = 10**n.polyval(ps,x)
n.savetxt("Proba-ELG-z-0.701838.dat", n.transpose([n.hstack((0,xb[sel][0]*0.9,xb[sel],xb[sel][-1]*1.1,1e4)), n.hstack((0,0,proba,0,0))]), header= "delta probability")

ps = n.loadtxt("fit-polynomial-ELG-z08.data")
proba = 10**n.polyval(ps,x)
n.savetxt("Proba-ELG-z-0.818843.dat", n.transpose([n.hstack((0,xb[sel][0]*0.9,xb[sel],xb[sel][-1]*1.1,1e4)), n.hstack((0,0,proba,0,0))]), header= "delta probability")

ps = n.loadtxt("fit-polynomial-LRG-z07.data")
proba = 10**n.polyval(ps,x)
n.savetxt("Proba-LRG-z-0.701838.dat", n.transpose([n.hstack((0,xb[sel][0]*0.9,xb[sel],xb[sel][-1]*1.1,1e4)), n.hstack((0,0,proba,0,0))]), header= "delta probability")

ps = n.loadtxt("fit-polynomial-LRG-z08.data")
proba = 10**n.polyval(ps,x)
n.savetxt("Proba-LRG-z-0.818843.dat", n.transpose([n.hstack((0,xb[sel][0]*0.9,xb[sel],xb[sel][-1]*1.1,1e4)), n.hstack((0,0,proba,0,0))]), header= "delta probability")


ps = n.loadtxt("fit-polynomial-QSO-z07.data")
proba = 10**n.polyval(ps,x)
n.savetxt("Proba-QSO-z-0.701838.dat", n.transpose([n.hstack((0,xb[sel][0]*0.9,xb[sel],xb[sel][-1]*1.1,1e4)), n.hstack((0,0,proba,0,0))]), header= "delta probability")

ps = n.loadtxt("fit-polynomial-QSO-z08.data")
proba = 10**n.polyval(ps,x)
n.savetxt("Proba-QSO-z-0.818843.dat", n.transpose([n.hstack((0,xb[sel][0]*0.9,xb[sel],xb[sel][-1]*1.1,1e4)), n.hstack((0,0,proba,0,0))]), header= "delta probability")


ps = n.loadtxt("fit-polynomial-QSO-z15.data")
proba = 10**n.polyval(ps,x)
n.savetxt("Proba-QSO-z-1.480160.dat", n.transpose([n.hstack((0,xb[sel][0]*0.9,xb[sel],xb[sel][-1]*1.1,1e4)), n.hstack((0,0,proba,0,0))]), header= "delta probability")
