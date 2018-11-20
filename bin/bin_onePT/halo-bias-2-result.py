import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
#import lib_functions_1pt as lib
import os
import sys

# data modules
import glob
import sys
import astropy.io.fits as fits
import os
from os.path import join
#import cPickle

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
delta_c = 1.686

mvir_dir = os.environ['MVIR_DIR']
#mvir_dir = join(os.environ['OBS_REPO'], 'data_2017/', 'mvir_dir')


import astropy.io.fits as fits
# loads summary file

from colossus.cosmology import cosmology
cosmology.setCosmology('planck15')

params = {'flat': True, 'H0': 67.77, 'Om0': 0.307115, 'Ob0': 0.048206, 'sigma8': 0.8228, 'ns': 0.9600, 'relspecies': False}
cosmology.setCosmology('myCosmo', params)


from colossus.lss import peaks
delta_c = peaks.collapseOverdensity(corrections = True, z = 0)

#delta_c = 1.686
nu = n.arange(0.2, 8.0, 0.01)
Masses = peaks.massFromPeakHeight(nu, 0.0)
sigma = ( delta_c / nu)



b_BH = lambda sigma, a, p, q:  1 + (a*(delta_c/sigma)**2. - q) / delta_c + (2*p/delta_c)/(1 + (a*(delta_c/sigma)**2.)**p)

qty = 'mvir'
#dir = join(os.environ['MD_DIR'])


s_low, s_high, s_mean, scale, b, bErr, volume, a, logmps = n.loadtxt(join(mvir_dir,  "halo-bias-measurement-summary.data"), unpack=True)

sel0 = (a==1) & (n.isnan(b)==False) 

bias = lambda sigma : b_BH(sigma, a=0.897, p=0.624, q=1.589)
alpha = lambda a : 0.346 - 0.059*a + 0.025*a**2.
logBeta = lambda a : 2.209 + 0.060 * a -0.021*a**2.

#m2sigma = interp1d(hmf.M, sigma )
#s_low = m2sigma( m_low )
#s_high = m2sigma( m_high )
#s_mean = m2sigma( m_mean )

bias_all = n.array([
b_BH(sigma, a=0.897+0.006, p=0.624+0.025, q=1.589+0.03), 
b_BH(sigma, a=0.897+0.006, p=0.624+0.025, q=1.589-0.03), 
b_BH(sigma, a=0.897+0.006, p=0.624-0.025, q=1.589-0.03), 
b_BH(sigma, a=0.897+0.006, p=0.624-0.025, q=1.589+0.03), 
b_BH(sigma, a=0.897-0.006, p=0.624-0.025, q=1.589-0.03), 
b_BH(sigma, a=0.897-0.006, p=0.624-0.025, q=1.589+0.03), 
b_BH(sigma, a=0.897-0.006, p=0.624+0.025, q=1.589-0.03),
b_BH(sigma, a=0.897-0.006, p=0.624+0.025, q=1.589+0.03)
])
bi_max = n.max(bias_all,axis = 0)
bi_min = n.min(bias_all,axis = 0)


bias_all = n.array([
b_BH(sigma, a=0.705+0.008, p=0.636+0.018, q=1.451+0.031), 
b_BH(sigma, a=0.705+0.008, p=0.636+0.018, q=1.451-0.031), 
b_BH(sigma, a=0.705+0.008, p=0.636-0.018, q=1.451-0.031), 
b_BH(sigma, a=0.705+0.008, p=0.636-0.018, q=1.451+0.031), 
b_BH(sigma, a=0.705-0.008, p=0.636-0.018, q=1.451-0.031), 
b_BH(sigma, a=0.705-0.008, p=0.636-0.018, q=1.451+0.031), 
b_BH(sigma, a=0.705-0.008, p=0.636+0.018, q=1.451-0.031),
b_BH(sigma, a=0.705-0.008, p=0.636+0.018, q=1.451+0.031)
])
bi2_max = n.max(bias_all,axis = 0)
bi2_min = n.min(bias_all,axis = 0)

x_data = n.log10(s_mean[sel0])
y_data = b[sel0]
y_data_err = bErr[sel0]

ps = n.array([0.897, 0.552, 1.589])
b_fun = lambda logsigma, a, p, q : b_BH(10**logsigma, a, p, q) 
print( "bias fit ----------------------"    )
pOpt, pCov=curve_fit(b_fun, x_data, y_data, ps, sigma = y_data_err, maxfev=50000000)#, bounds=boundaries)
chi2 = n.sum(((b_fun(x_data, pOpt[0], pOpt[1],pOpt[2])-y_data)/y_data_err)**2. ) 
ndof = (len(x_data) - len(ps)) 
print( "best params=", pOpt                                    )
print( "err=", pCov.diagonal()**0.5                            )
print( "chi2 ", chi2, ndof, chi2/ndof                          )
print( "P chi2 1-cdf", 1-stc2.cdf(int(chi2),ndof)              )
print( "---------------------------------------------------"   )
pOpt_ST01 = pOpt
pErr_ST01 = pCov.diagonal()**0.5

n.savetxt(join(mvir_dir,"biasFunction_parameters_ST01_MFonly_fit.txt"), n.transpose([pOpt_ST01, pErr_ST01]), header="A a p")


p.figure(2,(6,6))
p.axes([0.17, 0.17, 0.78, 0.78])

p.plot(-n.log10(sigma), b_BH(sigma, a=0.897, p=0.624, q=1.589),'k', lw=2)
p.fill_between(-n.log10(sigma), y1=bi_min, y2=bi_max, color='k',label='only HMF',alpha=0.2)

p.plot(-n.log10(sigma), b_fun(n.log10(sigma), pOpt[0], pOpt[1],pOpt[2]), 'm', lw=2)
p.fill_between(-n.log10(sigma), y1=bi2_min, y2=bi2_max, color='m',label='only bias',alpha=0.2)

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
p.xlim((-0.5,0.32))
gl=p.legend(loc=0, frameon=False)
p.grid()
p.savefig(join(mvir_dir,"halo_bias_v2.png"))
p.clf()
	