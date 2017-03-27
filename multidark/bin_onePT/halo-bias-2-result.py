import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
import lib_functions_1pt as lib
import os
import sys

# data modules
import glob
import sys
import astropy.io.fits as fits
import os
from os.path import join
import cPickle

# numerical modules
import numpy as n
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.interpolate import interp2d
from scipy.stats import norm
from scipy.interpolate import griddata
from scipy.stats import chi2 as stc2

# plotting modules
import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p

# mass function theory
from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)


qty = 'mvir'
dir = join(os.environ['MD_DIR'])
s_low, s_high, s_mean, scale, b, bErr, volume, a, logmps = n.loadtxt(join(os.environ['MVIR_DIR'],  "halo-bias-measurement-summary.data"), unpack=True)

sel0 = (a==1) & (n.isnan(b)==False) 

bias = lambda sigma : lib.b_BH(sigma, a=0.8915, p=0.5524, q=1.578)
alpha = lambda a : 0.346 - 0.059*a + 0.025*a**2.
logBeta = lambda a : 2.209 + 0.060 * a -0.021*a**2.

#m2sigma = interp1d(lib.hmf.M, lib.hmf.sigma )
#s_low = m2sigma( m_low )
#s_high = m2sigma( m_high )
#s_mean = m2sigma( m_mean )

bias_all = n.array([
lib.b_BH(lib.hmf.sigma, a=0.8915+0.01, p=0.5524+0.05, q=1.578+0.07), 
lib.b_BH(lib.hmf.sigma, a=0.8915+0.01, p=0.5524+0.05, q=1.578-0.07), 
lib.b_BH(lib.hmf.sigma, a=0.8915+0.01, p=0.5524-0.05, q=1.578-0.07), 
lib.b_BH(lib.hmf.sigma, a=0.8915+0.01, p=0.5524-0.05, q=1.578+0.07), 
lib.b_BH(lib.hmf.sigma, a=0.8915-0.01, p=0.5524-0.05, q=1.578-0.07), 
lib.b_BH(lib.hmf.sigma, a=0.8915-0.01, p=0.5524-0.05, q=1.578+0.07), 
lib.b_BH(lib.hmf.sigma, a=0.8915-0.01, p=0.5524+0.05, q=1.578-0.07),
lib.b_BH(lib.hmf.sigma, a=0.8915-0.01, p=0.5524+0.05, q=1.578+0.07)
])
bi_max = n.max(bias_all,axis = 0)
bi_min = n.min(bias_all,axis = 0)


bias_all = n.array([
lib.b_BH(lib.hmf.sigma, a=0.7400+0.0084, p=0.6140+0.0184, q=1.6468+0.0259), 
lib.b_BH(lib.hmf.sigma, a=0.7400+0.0084, p=0.6140+0.0184, q=1.6468-0.0259), 
lib.b_BH(lib.hmf.sigma, a=0.7400+0.0084, p=0.6140-0.0184, q=1.6468-0.0259), 
lib.b_BH(lib.hmf.sigma, a=0.7400+0.0084, p=0.6140-0.0184, q=1.6468+0.0259), 
lib.b_BH(lib.hmf.sigma, a=0.7400-0.0084, p=0.6140-0.0184, q=1.6468-0.0259), 
lib.b_BH(lib.hmf.sigma, a=0.7400-0.0084, p=0.6140-0.0184, q=1.6468+0.0259), 
lib.b_BH(lib.hmf.sigma, a=0.7400-0.0084, p=0.6140+0.0184, q=1.6468-0.0259),
lib.b_BH(lib.hmf.sigma, a=0.7400-0.0084, p=0.6140+0.0184, q=1.6468+0.0259)
])
bi2_max = n.max(bias_all,axis = 0)
bi2_min = n.min(bias_all,axis = 0)

x_data = n.log10(s_mean[sel0])
y_data = b[sel0]
y_data_err = bErr[sel0]

ps = n.array([0.8915, 0.552, 1.578])
b_fun = lambda logsigma, a, p, q : lib.b_BH(10**logsigma, a, p, q) 
print "bias fit ----------------------"
pOpt, pCov=curve_fit(b_fun, x_data, y_data, ps, sigma = y_data_err, maxfev=50000000)#, bounds=boundaries)
chi2 = n.sum(((b_fun(x_data, pOpt[0], pOpt[1],pOpt[2])-y_data)/y_data_err)**2. ) 
ndof = (len(x_data) - len(ps)) 
print "best params=", pOpt
print "err=", pCov.diagonal()**0.5
print "chi2 ", chi2, ndof, chi2/ndof
print "P chi2 1-cdf", 1-stc2.cdf(int(chi2),ndof)
print "---------------------------------------------------"
pOpt_ST01 = pOpt
pErr_ST01 = pCov.diagonal()**0.5

n.savetxt(join(os.environ['MVIR_DIR'],"biasFunction_parameters_ST01_MFonly_fit.txt"), n.transpose([pOpt_ST01, pErr_ST01]), header="A a p")


p.figure(2,(6,6))
p.axes([0.17, 0.17, 0.78, 0.78])

p.plot(-n.log10(lib.hmf.sigma), lib.b_BH(lib.hmf.sigma, a=0.8915, p=0.5524, q=1.578),'k', lw=2)
p.fill_between(-n.log10(lib.hmf.sigma), y1=bi_min, y2=bi_max, color='k',label='only HMF',alpha=0.2)

p.plot(-n.log10(lib.hmf.sigma), b_fun(n.log10(lib.hmf.sigma), pOpt[0], pOpt[1],pOpt[2]), 'm', lw=2)
p.fill_between(-n.log10(lib.hmf.sigma), y1=bi2_min, y2=bi2_max, color='m',label='only bias',alpha=0.2)

sel = (sel0) & (volume==400.**3.)
p.errorbar(-n.log10(s_mean[sel]), b[sel], xerr=abs(n.log10(s_high[sel]) -n.log10(s_low[sel]) )/2., yerr = bErr[sel], rasterized=True, fmt='none', label='M04')

sel = (sel0) & (volume==1000.**3.)
p.errorbar(-n.log10(s_mean[sel]), b[sel], xerr=abs(n.log10(s_high[sel]) -n.log10(s_low[sel]) )/2., yerr = bErr[sel], rasterized=True, fmt='none', label='M10')

sel = (sel0) & (volume==2500.**3.)
p.errorbar(-n.log10(s_mean[sel]), b[sel], xerr=abs(n.log10(s_high[sel]) -n.log10(s_low[sel]) )/2., yerr = bErr[sel], rasterized=True, fmt='none', label='M25 or M25n')

sel = (sel0) & (volume==4000.**3.)
p.errorbar(-n.log10(s_mean[sel]), b[sel], xerr=abs(n.log10(s_high[sel]) -n.log10(s_low[sel]) )/2., yerr = bErr[sel], rasterized=True, fmt='none', label='M40 or M40n')

p.xlabel(r'$\log_{10}(\sigma^{-1})$')
p.ylabel(r'bias')
p.ylim((0.5,5))
p.xlim((-0.5,0.3))
gl=p.legend(loc=0, frameon=False)
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],"halo_bias_v2.png"))
p.clf()
	