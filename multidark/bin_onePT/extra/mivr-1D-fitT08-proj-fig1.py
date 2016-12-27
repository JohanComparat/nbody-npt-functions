import glob
import sys
import cPickle
from os.path import join
import numpy as n
import astropy.io.fits as fits
import os

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu

import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

from scipy.optimize import minimize
from scipy.optimize import curve_fit

from scipy.interpolate import interp1d
from scipy.misc import derivative

delta_c = 1.686

f_SMT = lambda sigma, A, a, p: A * (2.*a/n.pi)**(0.5) * ( 1 + ((delta_c/sigma)**2./a) **(p) ) * n.e**( - a * (delta_c/sigma)**2. / 2.) * (delta_c/sigma)

log_f_ST01 = lambda logSigma, A, a, p : n.log10( f_SMT(10.**logSigma, A, a, p) )

outfile=open(join("..","mvir","mvir-cen-diff-function-z0-params.pkl"), 'r')
outCF = cPickle.load(outfile)
outfile.close()
pOpt, pCov = outCF
A0, a0, p0 = outCF[0]

rhom = cosmo.critical_density(0.).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(0.)/(100*uu.km/(uu.Mpc*uu.s)))**1.

frho = lambda sigma : f_SMT(sigma, A0, a0, p0) * rhom.value

X = n.arange(-0.4, 0.6, 0.01)

msigmaFile=join("..", "Pk_DM_CLASS", "hmf_highz_medz_lowz_planck", "mVector_z_0.0.txt")
DATA = n.loadtxt(msigmaFile,unpack=True)

M=DATA[0]
sigma = DATA[1]
m2sigma = interp1d(M, sigma)
toderive = interp1d(n.log(M), DATA[2])
dlnsigdm = interp1d(M[10:-10], derivative(toderive, n.log(M[10:-10])))

mass = M[10:-10]
jac =  dlnsigdm(mass) 
sig = m2sigma(mass)
counts = frho(sig) * jac / mass * 1e9

volumes = 10**n.arange(5,13.1,0.1)
mass100=n.ones_like(volumes)*-1
mass10k=n.ones_like(volumes)*-1
mass1M=n.ones_like(volumes)*-1
for ii, volume in enumerate(volumes): # = 1e9
	counts = frho(sig) * jac / mass * volume
	counts_c = n.array([ n.sum(counts[jj:]) for jj in range(len(counts)) ]) 
	c2m = interp1d(counts_c, n.log10(mass))
	print ii, n.min(counts_c), n.max(counts_c)
	if n.min(counts_c)<100 :
		mass100[ii] = c2m(100)
		mass10k[ii] = c2m(10000)
		mass1M[ii] = c2m(1000000)
	if n.min(counts_c)<10000 :
		mass10k[ii] = c2m(10000)
		mass1M[ii] = c2m(1000000)
	if n.min(counts_c)<1000000 :
		mass1M[ii] = c2m(1000000)

sys.exit()
		
n.savetxt(join("..",qty,"volume-number.txt"), n.transpose([n.arange(5,13.1,0.1), mass100, mass10k, mass1M]), fmt='%2.2f', header=" logVol massN100 massN10k massN1M ")


sys.exit()

# redshift selection
zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
# mass selection
mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0])) &(data["log_"+qty+"_max"]<n.log10(limits_04[1]))) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])) &(data["log_"+qty+"_max"]<n.log10(limits_10[1]))) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>n.log10(limits_25[0])) &(data["log_"+qty+"_max"]<n.log10(limits_25[1]))) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>n.log10(limits_40[0]))&(data["log_"+qty+"_max"]<n.log10(limits_40[1]))) 

# minimum number counts selection
nSel = (data['dN_counts_'+cos]>NminCount)
# altogether
ok = (zSel) & (mSel) & (nSel)

# x coordinates
logsigM1 = n.log10(1./data['sigmaM'][ok])#
sig = data['sigmaM'][ok]
#print n.min(logsigM1), n.max(logsigM1)
log_m200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
m200c = 10**log_m200c
# mean density array normalization
rhom = cosmo.critical_density(data["redshift"][ok]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][ok])/(100*uu.km/(uu.Mpc*uu.s)))**1.

x_data = n.log10(data['sigmaM'][ok])#
z_data = data['redshift'][ok]#
f_data = (m200c * data["dNdVdlnM_"+cos][ok]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"][ok]) )

jac =  abs(data["dlnsigmaM1_o_dlnM"][ok])
ids = n.argsort(m200c)
sigma = sig[ids]
jacobian = jac[ids]
mass = m200c[ids]

volumes = 10**n.arange(5,13.1,0.1)
mass100=n.ones_like(volumes)*-1
mass1000=n.ones_like(volumes)*-1
mass10000=n.ones_like(volumes)*-1
mass100k=n.ones_like(volumes)*-1
for ii, volume in enumerate(volumes): # = 1e9
	counts = frho(sigma) * jacobian / mass * volume
	counts_c = n.array([ n.sum(counts[jj:]) for jj in range(len(counts)) ]) 
	c2m = interp1d(counts_c, n.log10(mass))
	print ii, n.min(counts_c), n.max(counts_c)
	if n.min(counts_c)<100 :
		mass100[ii] = c2m(100)
		mass1000[ii] = c2m(1000)
		mass10000[ii] = c2m(10000)
		mass100k[ii] = c2m(100000)
	if n.min(counts_c)<1000 :
		mass1000[ii] = c2m(1000)
		mass10000[ii] = c2m(10000)
		mass100k[ii] = c2m(100000)
	if n.min(counts_c)<10000 :
		mass10000[ii] = c2m(10000)
		mass100k[ii] = c2m(100000)
	if n.min(counts_c)<100000 :
		mass100k[ii] = c2m(100000)

n.savetxt(join("..",qty,"volume-number.txt"), n.transpose([n.arange(5,13.1,0.1), mass100,  mass1000,  mass10000,  mass100k]), fmt='%2.2f', header=" logVol massN100 massN1k massN10k massN100k ")
