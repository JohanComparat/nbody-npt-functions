import glob
import sys
import astropy.io.fits as fits
import os
from os.path import join
import cPickle
import time
# numerical modules
import numpy as n
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import curve_fit

# plotting modules
import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

# mass function theory
from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

from scipy.interpolate import interp1d
from scipy.integrate import quad

# Fitting functions
# velocity function
vf = lambda v, A, v0, alpha, beta : n.log10( 10**A * (10**v/10**v0)**(-beta) * n.e**(- (10**v/10**v0)**(alpha) ) )
"""
# sheth and tormen function
fnu_ST = lambda nu, A, a, p: A * (2./n.pi)**(0.5) * (a*nu) * ( 1 + (a * nu) **(-2.*p) ) * n.e**( -( a * nu )**2. / 2.) 
log_fnu_ST = lambda logNu, Anorm, a, q : n.log10( fnu_ST(10.**logNu, Anorm, a, q) )
log_fnu_ST_ps = lambda logNu, ps : n.log10( fnu_ST(10.**logNu, ps[0], ps[1], ps[2]) )

p_ST_despali = [0.333, 0.794, 0.247]
p_ST_sheth = [0.3222, 0.707, 0.3]
p_T08_klypin = [0.224, 1.67, 1.80, 1.48]
p_T08_RP = [0.144, 1.351, 3.113, 1.187]


# tinker function 
fsigma_T08 = lambda sigma, A, a, b, c :  A *(1+ (b/sigma)**a) * n.e**(-c/sigma**2.)
log_fsigma_T08 = lambda  logsigma, A, a, b, c :  n.log10(A *(1+ (b/(10**logsigma))**a) * n.e**(-c/(10**logsigma)**2.))
log_fsigma_T08_ps = lambda logsigma, ps : log_fsigma_T08(logsigma, ps[0], ps[1], ps[2], ps[3])

#fnu_SMT = lambda nu, Anorm, a, q : Anorm * (2./n.pi)**(0.5) * (a*nu) * ( 1 + (a*nu) **(-2*q) ) * n.e**( - (a*nu)**2. / 2.) # 
#fnu_SMT = lambda nu, Anorm, a, q : Anorm *a * (2.*n.pi)**(-0.5) *  ( 1 + (a**2*nu) **(-q) ) * n.e**( - a**2*nu / 2.)
#log_fnu_ST01 = lambda logNu, Anorm, a, q : n.log10( fnu_SMT(10.**logNu, Anorm, a, q) )
#log_fnu_ST01_ps = lambda logNu, ps : n.log10( fnu_SMT(10.**logNu, ps[0], ps[1], ps[2]) )
"""


delta_c = 1.686
sigma = n.arange(0.05,10,0.05)
X = n.arange(-0.6, 0.5, 0.01) #n.log10(1./sigma)
sigma = 10**-X 

f_dsp_nu = lambda nu, A, a, p: A* (a / (nu * 2 * n.pi))**(0.5) * ( 1 + 1. / (a*nu) **p ) * n.e**( - a * nu / 2.)
nufnu_dsp = lambda nu, A, a, p: A* ((2 * a * nu) / ( n.pi))**(0.5) * ( 1 + (a*nu) **(-p) ) * n.e**( - a * nu / 2.)

# nu = (delta_c/sigma)**2.

#f_dsp = lambda sigma, A, a, p : A* (a / ((delta_c/sigma)**2. * 2 * n.pi))**(0.5) * ( 1 + 1. / (a*(delta_c/sigma)**2.) **p ) * n.e**( - a * (delta_c/sigma)**2. / 2.)
f_dsp = lambda sigma, A, a, p: A * ((2. * a * (delta_c/sigma)**2.) / (  n.pi))**(0.5) * ( 1 + (a*(delta_c/sigma)**2.) **(-p) ) * n.e**( - a * (delta_c/sigma)**2. / 2.)

f_T08 = lambda sigma, A, a, b, c : A*( (sigma/b)**(-a) + 1 )*n.e**(-c/sigma**2.)
f_ST = lambda sigma, A, a, p: A* (2*a/n.pi)**0.5 * ( 1 + (delta_c**2./(a*sigma**2.))**(p) )*(delta_c/sigma)*n.e**(-a*delta_c**2./(2.*sigma**2.))
f_BH = lambda sigma, A, a, p, q: A* (2./n.pi)**0.5 * ( 1 + (sigma**2./(a**delta_c*2.))**(p) )*(delta_c*a**0.5/sigma)**(q)*n.e**(-a*delta_c**2./(2.*sigma**2.))
b_BH = lambda sigma, a, p, q:  1 + (a*(delta_c/sigma)**2. - q) / delta_c + (2*p/delta_c)/(1 + (a*(delta_c/sigma)**2.))**p


"""
M200c
ftT08 = f_T08(sigma, 0.186, 1.47, 2.57, 1.19) 
ftSk14 = f_T08(sigma, 0.18587, 1.46690, 2.57110, 1.19396)
ftK16 = f_T08(sigma, 0.224, 1.67, 1.80, 1.48) 
ftA12 = f_T08 (sigma, 0.201, 1.7, 2.08, 1.172)
ftW13 = f_T08 (sigma, 0.282, 2.163, 1.406, 1.21)
ftST01 = f_ST(sigma, 0.3222, 0.707, 0.3)
#ftD16 = f_dsp(sigma, 0.287, 0.903, 0.322 )
ftC16 = f_T08(sigma, 0.12, 1.19, 3.98, 1.35)
ftC16st = f_dsp(sigma, 0.2906, 0.8962, 0.1935 )
"""
ftC16 = f_BH(sigma, 0.279, 0.908, 0.671, 1.737)
ftC16st = f_dsp(sigma, 0.317, 0.819, 0.113)

ftD16 = f_dsp(sigma, 0.333, 0.794, 0.247 )
ftRP16 = f_T08 (sigma, 0.144, 1.351, 3.113, 1.187)
ftT08 = f_T08(sigma, 0.200, 1.47, 2.57, 1.19) # 300 at z=0
#ftST01 = f_ST(sigma, 0.3222, 0.707, 0.3)
ftST02 = f_ST(sigma, 0.3222, 0.75, 0.3)
ftBH11 = f_BH(sigma, 0.333, 0.788, 0.807, 1.795)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.plot(X, ftC16/ftD16, label='this work (14)', lw=1.5)
p.plot(X, ftC16st/ftD16, label='this work (15)', lw=1.5)
#p.plot(X, ftSk14/ftT08, label='Skillman 14', lw=1.5)
p.plot(X, ftRP16/ftD16, label='Rodriguez-Puebla 16', lw=1.5)
#p.plot(X, ftK16/ftD16, label='Klypin 16', lw=1.5)
#p.plot(X, ftA12/ftT08, label='Angulo 12', ls='--')
#p.plot(X, ftW13/ftT08, label='Watson 13', ls='--')
p.plot(X, ftST02/ftD16, label='Sheth Tormen 02', lw=1.5)
p.plot(X, ftT08/ftD16, label='Tinker 08', lw=1.5)
p.plot(X, ftBH11/ftD16, label='Bhattacharya 11', lw=1.5)
p.axhline(1.10, ls='-.', color='k', label='10%')
p.axhline(0.9, ls='-.', color='k')
p.ylabel('model / Despali 2016')
p.xlabel(r'$log(1/\sigma)$')
gl = p.legend(loc=0,fontsize=12)
gl.set_frame_on(False)
p.title(r'$f(\sigma)$ comparison M$_{vir}$')
p.grid()
p.ylim((0,2))
p.xlim((-0.6,0.4))
p.savefig(join(os.environ['MVIR_DIR'],"mvir-fit-comparison-D16.png"))
p.clf()

sys.exit()

p.plot(lgsig2m(X), ftC16st/ftT08, label='c16 st01', lw=1.5)
p.plot(lgsig2m(X), ftST01/ftT08, label='Sheth 01', ls='--')
p.plot(lgsig2m(X), ftT08/ftT08, label='Tinker 08', ls='--')
#p.plot(lgsig2m(X), ftA12/ftT08, label='Angulo 12', ls='--')
#p.plot(lgsig2m(X), ftW13/ftT08, label='Watson 13', ls='--')
#p.plot(lgsig2m(X), ftSk14/ftT08, label='Skillman 14', lw=1.5)
p.plot(lgsig2m(X), ftD16/ftT08, label='Despali 16', lw=1.5)
p.plot(lgsig2m(X), ftK16/ftT08, label='Klypin 16', lw=1.5)
p.plot(lgsig2m(X), ftRP16/ftT08, label='Rodriguez-Puebla 16', lw=1.5)
p.plot(lgsig2m(X), ftC16/ftT08, label='c16 t08', lw=1.5)
p.axhline(1.10, ls='-.', color='k', label='10%')
p.axhline(0.9, ls='-.', color='k')
p.ylabel(r'$f(\sigma)$ / model(Tinker 2008)')
p.xlabel(r'$log_{10}(M_{200c})$')
gl = p.legend(loc=0,fontsize=12)
gl.set_frame_on(False)
p.grid()
p.ylim((0,2))
p.xlim(( lgsig2m(-0.6), lgsig2m(0.4) ))
p.savefig(join(os.environ['MVIR_DIR'],"mvir-fit-comparison-T08-xmass.png"))
p.clf()


p.plot(lgsig2m(X), ftC16/ftK16, label='c16 t08', lw=1.5)
p.plot(lgsig2m(X), ftC16st/ftK16, label='c16 st01', lw=1.5)
p.plot(lgsig2m(X), ftST01/ftK16, label='Sheth 01', ls='--')
p.plot(lgsig2m(X), ftT08/ftK16, label='Tinker 08', ls='--')
#p.plot(lgsig2m(X), ftA12/ftK16, label='Angulo 12', ls='--')
#p.plot(lgsig2m(X), ftW13/ftK16, label='Watson 13', ls='--')
#p.plot(lgsig2m(X), ftSk14/ftK16, label='Skillman 14', lw=1.5)
p.plot(lgsig2m(X), ftK16/ftK16, label='Klypin 16', lw=1.5)
p.plot(lgsig2m(X), ftD16/ftK16, label='Despali 16', lw=1.5)
p.plot(lgsig2m(X), ftRP16/ftK16, label='Rodriguez-Puebla 16', lw=1.5)
p.axhline(1.10, ls='-.', color='k', label='10%')
p.axhline(0.9, ls='-.', color='k')
p.ylabel(r'$f(\sigma)$ / model(Klypin 2016)')
p.xlabel(r'$log_{10}(M_{200c})$')
gl = p.legend(loc=0,fontsize=12)
gl.set_frame_on(False)
#p.title(r'$f(\sigma)$ comparison M$_{vir}$ Planck 2014 Cosmology')
p.grid()
p.ylim((0,2))
p.xlim(( lgsig2m(-0.6), lgsig2m(0.4) ))
p.savefig(join(os.environ['MVIR_DIR'],"mvir-fit-comparison-K16-xMass.png"))
p.clf()


p.plot(X, ftC16/ftK16, label='c16 t08', lw=1.5)
p.plot(X, ftC16st/ftK16, label='c16 st01', lw=1.5)
p.plot(X, ftRP16/ftK16, label='Rodriguez-Puebla 16', lw=1.5)
p.plot(X, ftD16/ftK16, label='Despali 16', lw=1.5)
p.axhline(1.10, ls='-.', color='k', label='10%')
p.axhline(0.9, ls='-.', color='k')

p.ylabel('model / Klypin 2016')
p.xlabel(r'$log(1/\sigma)$')
gl = p.legend(loc=0,fontsize=12)
gl.set_frame_on(False)
p.title(r'$f(\sigma)$ comparison M$_{vir}$ Planck 2014 Cosmology')
p.grid()
p.ylim((0,2))
p.xlim((-0.6,0.4))
p.savefig(join(os.environ['MVIR_DIR'],"mvir-fit-comparison-K16.png"))
p.clf()