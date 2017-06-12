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
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p

from scipy.optimize import minimize
from scipy.optimize import curve_fit

from scipy.interpolate import interp1d
from scipy.misc import derivative

msigmaFile=join("..", "Pk_DM_CLASS", "hmf_highz_medz_lowz_planck", "mVector_z_0.0.txt")
DATA = n.loadtxt(msigmaFile,unpack=True)
lgsig2m = interp1d(n.log(1/DATA[1]), n.log10(DATA[0]))
ftM13 = interp1d(DATA[1], DATA[4] )

delta_c = 1.686
sigma = n.arange(0.05,10,0.05)
X = n.arange(-0.6, 0.4, 0.01) #n.log10(1./sigma)
sigma = 10**-X 

f_dsp_nu = lambda nu, A, a, p: A* (a / (nu * 2 * n.pi))**(0.5) * ( 1 + 1. / (a*nu) **p ) * n.e**( - a * nu / 2.)
nufnu_dsp = lambda nu, A, a, p: A* ((2 * a * nu) / ( n.pi))**(0.5) * ( 1 + (a*nu) **(-p) ) * n.e**( - a * nu / 2.)

# nu = (delta_c/sigma)**2.

#f_dsp = lambda sigma, A, a, p : A* (a / ((delta_c/sigma)**2. * 2 * n.pi))**(0.5) * ( 1 + 1. / (a*(delta_c/sigma)**2.) **p ) * n.e**( - a * (delta_c/sigma)**2. / 2.)
f_dsp = lambda sigma, A, a, p: A * ((2. * a * (delta_c/sigma)**2.) / (  n.pi))**(0.5) * ( 1 + (a*(delta_c/sigma)**2.) **(-p) ) * n.e**( - a * (delta_c/sigma)**2. / 2.)

f_T08 = lambda sigma, A, a, b, c : A*( (sigma/b)**(-a) + 1 )*n.e**(-c/sigma**2.)
f_ST = lambda sigma, A, a, p: A* (2*a/n.pi)**0.5 * ( 1 + (delta_c**2./(a*sigma**2.))**(p) )*(delta_c/sigma)*n.e**(-a*delta_c**2./sigma**2./2.)

ftT08 = f_T08(sigma, 0.186, 1.47, 2.57, 1.19) 
ftSk14 = f_T08(sigma, 0.18587, 1.46690, 2.57110, 1.19396)
ftC16 = f_T08(sigma, 0.12, 1.19, 3.98, 1.35)
ftK16 = f_T08(sigma, 0.224, 1.67, 1.80, 1.48) 
ftRP16 = f_T08 (sigma, 0.144, 1.351, 3.113, 1.187)
ftA12 = f_T08 (sigma, 0.201, 1.7, 2.08, 1.172)
ftW13 = f_T08 (sigma, 0.282, 2.163, 1.406, 1.21)
ftST01 = f_ST(sigma, 0.3222, 0.707, 0.3)
ftD16 = f_dsp(sigma, 0.287, 0.903, 0.322 )
ftC16st = f_dsp(sigma, 0.2906, 0.8962, 0.1935 )

p.plot(X, ftC16/ftT08, label='c16 t08', lw=2)
p.plot(X, ftC16st/ftT08, label='c16 st01', lw=2)
#p.plot(X, ftSk14/ftT08, label='Skillman 14', lw=2)
p.plot(X, ftRP16/ftT08, label='Rodriguez-Puebla 16', lw=2)
p.plot(X, ftK16/ftT08, label='Klypin 16', lw=2)
#p.plot(X, ftA12/ftT08, label='Angulo 12', ls='--')
#p.plot(X, ftW13/ftT08, label='Watson 13', ls='--')
p.plot(X, ftST01/ftT08, label='Sheth Tormen 01', ls='--')
p.plot(X, ftD16/ftT08, label='Despali 16', lw=2)
p.axhline(1.10, ls='-.', color='k', label='10%')
p.axhline(0.9, ls='-.', color='k')
p.ylabel('model / Tinker 2008')
p.xlabel(r'$log(1/\sigma)$')
gl = p.legend(loc=0,fontsize=12)
gl.set_frame_on(False)
p.title(r'$f(\sigma)$ comparison M$_{200c}$')
p.grid()
p.ylim((0,2))
p.xlim((-0.6,0.4))
p.savefig(join("..","M200c","M200c-fit-comparison-T08.png"))
p.clf()

p.plot(lgsig2m(X), ftC16st/ftT08, label='c16 st01', lw=2)
p.plot(lgsig2m(X), ftST01/ftT08, label='Sheth 01', ls='--')
p.plot(lgsig2m(X), ftT08/ftT08, label='Tinker 08', ls='--')
#p.plot(lgsig2m(X), ftA12/ftT08, label='Angulo 12', ls='--')
#p.plot(lgsig2m(X), ftW13/ftT08, label='Watson 13', ls='--')
#p.plot(lgsig2m(X), ftSk14/ftT08, label='Skillman 14', lw=2)
p.plot(lgsig2m(X), ftD16/ftT08, label='Despali 16', lw=2)
p.plot(lgsig2m(X), ftK16/ftT08, label='Klypin 16', lw=2)
p.plot(lgsig2m(X), ftRP16/ftT08, label='Rodriguez-Puebla 16', lw=2)
p.plot(lgsig2m(X), ftC16/ftT08, label='c16 t08', lw=2)
p.axhline(1.10, ls='-.', color='k', label='10%')
p.axhline(0.9, ls='-.', color='k')
p.ylabel(r'$f(\sigma)$ / model(Tinker 2008)')
p.xlabel(r'$log_{10}(M_{200c})$')
gl = p.legend(loc=0,fontsize=12)
gl.set_frame_on(False)
p.grid()
p.ylim((0,2))
p.xlim(( lgsig2m(-0.6), lgsig2m(0.4) ))
p.savefig(join("..","M200c","M200c-fit-comparison-T08-xmass.png"))
p.clf()


p.plot(lgsig2m(X), ftC16/ftK16, label='c16 t08', lw=2)
p.plot(lgsig2m(X), ftC16st/ftK16, label='c16 st01', lw=2)
p.plot(lgsig2m(X), ftST01/ftK16, label='Sheth 01', ls='--')
p.plot(lgsig2m(X), ftT08/ftK16, label='Tinker 08', ls='--')
#p.plot(lgsig2m(X), ftA12/ftK16, label='Angulo 12', ls='--')
#p.plot(lgsig2m(X), ftW13/ftK16, label='Watson 13', ls='--')
#p.plot(lgsig2m(X), ftSk14/ftK16, label='Skillman 14', lw=2)
p.plot(lgsig2m(X), ftK16/ftK16, label='Klypin 16', lw=2)
p.plot(lgsig2m(X), ftD16/ftK16, label='Despali 16', lw=2)
p.plot(lgsig2m(X), ftRP16/ftK16, label='Rodriguez-Puebla 16', lw=2)
p.axhline(1.10, ls='-.', color='k', label='10%')
p.axhline(0.9, ls='-.', color='k')
p.ylabel(r'$f(\sigma)$ / model(Klypin 2016)')
p.xlabel(r'$log_{10}(M_{200c})$')
gl = p.legend(loc=0,fontsize=12)
gl.set_frame_on(False)
#p.title(r'$f(\sigma)$ comparison M$_{200c}$ Planck 2014 Cosmology')
p.grid()
p.ylim((0,2))
p.xlim(( lgsig2m(-0.6), lgsig2m(0.4) ))
p.savefig(join("..","M200c","M200c-fit-comparison-K16-xMass.png"))
p.clf()


p.plot(X, ftC16/ftK16, label='c16 t08', lw=2)
p.plot(X, ftC16st/ftK16, label='c16 st01', lw=2)
p.plot(X, ftRP16/ftK16, label='Rodriguez-Puebla 16', lw=2)
p.plot(X, ftD16/ftK16, label='Despali 16', lw=2)
p.axhline(1.10, ls='-.', color='k', label='10%')
p.axhline(0.9, ls='-.', color='k')

p.ylabel('model / Klypin 2016')
p.xlabel(r'$log(1/\sigma)$')
gl = p.legend(loc=0,fontsize=12)
gl.set_frame_on(False)
p.title(r'$f(\sigma)$ comparison M$_{200c}$ Planck 2014 Cosmology')
p.grid()
p.ylim((0,2))
p.xlim((-0.6,0.4))
p.savefig(join("..","M200c","M200c-fit-comparison-K16.png"))
p.clf()