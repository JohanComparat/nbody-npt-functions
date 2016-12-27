from os.path import join
import numpy as n
import astropy.io.fits as fits
import sys
import os
import lib_functions_1pt as lib
import cPickle

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu

import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

from scipy.interpolate import interp1d
from scipy.misc import derivative
import astropy.units as uu

from scipy.optimize import curve_fit
#Quantity studied
qty = "mvir"
# working directory
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits"))[1].data

NminCount = 1000
logNpmin = 3

zmin = -0.01
zmax = 0.001


# x coordinates definition
logsig = -n.log10(data['sigmaM'])#
lognu = n.log10(data['nu2']**0.5)
#log_mvir = data["log_"+qty]
log_mvir = data["log_"+qty] - n.log10(cosmo.h)
mvir = 10**data["log_"+qty] / cosmo.h

#=======================
#=======================
cos = 'cen'
#=======================
#=======================
# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )
# error on y position
#=================
error = data["dN_counts_"+cos]**(-0.5)+0.005


# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, logNpmin)
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelCen)
# selection per box :
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc_new_rockS')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')


p0 = [0.328, 0.8, 0.2]
sigma = n.arange(0.05,10,0.05)

"""
p.plot(data['nu2'][ok], log_MF[ok], 'bo')
p.plot(data['nu2'][ok], log_fnu_ST_ps(n.log10(data['nu2'][ok]**0.5), p0), 'k+')
p.xscale('log')
p.ylim((-4.5,0))
p.grid()
p.savefig(join(dir,"fit-"+cos+"-first-guess.png"))
p.clf()
"""

x_data = n.log10(data['nu2'][ok]**0.5)
y_data = log_MF[ok]
y_err = error[ok]

pOptCen, pCovCen = lib.fit_mvir_function_z0(data[ok], x_data , y_data , y_err , p0, tolerance = 0.03, cos = cos)


#=======================
#=======================
cos = 'sat'
#=======================
#=======================
# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )
# error on y position
#=================
error = data["dN_counts_"+cos]**(-0.5)+0.005


# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, logNpmin)
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelCen)
# selection per box :
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc_new_rockS')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')


p0 = [0.328, 0.8, 0.2]
sigma = n.arange(0.05,10,0.05)

"""
p.plot(data['nu2'][ok], log_MF[ok], 'bo')
p.plot(data['nu2'][ok], log_fnu_ST_ps(n.log10(data['nu2'][ok]**0.5), p0), 'k+')
p.xscale('log')
p.ylim((-4.5,0))
p.grid()
p.savefig(join(dir,"fit-"+cos+"-first-guess.png"))
p.clf()
"""

x_data = n.log10(data['nu2'][ok]**0.5)
y_data = log_MF[ok]
y_err = error[ok]

pOptSat, pCovSat = lib.fit_mvir_function_z0(data[ok], x_data , y_data , y_err , p0, tolerance = 0.03, cos = cos)

# deviations from universality centrals

zmin = -0.01
zmax = 2.5

cos = 'sat'
# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )
# error on y position
#=================
error = data["dN_counts_"+cos]**(-0.5)+0.005


# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, logNpmin)
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelCen)


x_data = n.log10(data['nu2'][ok]**0.5)
y_data = log_MF[ok]
y_err = error[ok]

y_model=lib.log_fnu_ST_ps(x_data, pOptSat)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.errorbar(x_data, 10**(y_data-y_model), yerr = y_err , rasterized=True, fmt='none')

p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
p.xlabel(r'$log_{10}(\nu)$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
#p.xlim((-0.7,0.6))
p.ylim((0.5,1.5))
#p.yscale('log')
p.grid()
p.savefig(join(dir,"fit-"+cos+"-differential-function-universal-residual-log.png"))
p.clf()

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(x_data, 10**(y_data-y_model), c=data['redshift'][ok], s=5, marker='o',rasterized=True)#, vmin=zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")

p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
p.xlabel(r'$log_{10}(\nu)$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
#p.xlim((-0.7,0.6))
p.ylim((0.5,1.5))
#p.yscale('log')
p.grid()
p.savefig(join(dir,"fit-"+cos+"-differential-function-universal-residual-log-cz.png"))
p.clf()
# deviations from universality centrals

zmin = -0.01
zmax = 2.5

cos = 'cen'
# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )
# error on y position
#=================
error = data["dN_counts_"+cos]**(-0.5)+0.005


# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, logNpmin)
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelCen)


x_data = n.log10(data['nu2'][ok]**0.5)
y_data = log_MF[ok]
y_err = error[ok]

y_model=lib.log_fnu_ST_ps(x_data, pOptCen)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.errorbar(x_data, 10**(y_data-y_model), yerr = y_err , rasterized=True, fmt='none')

p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
p.xlabel(r'$log_{10}(\nu)$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
#p.xlim((-0.7,0.6))
p.ylim((0.5,1.5))
#p.yscale('log')
p.grid()
p.savefig(join(dir,"fit-"+cos+"-differential-function-universal-residual-log.png"))
p.clf()

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(x_data, 10**(y_data-y_model), c=data['redshift'][ok], s=5, marker='o',rasterized=True)#, vmin=zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")

p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
p.xlabel(r'$log_{10}(\nu)$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
#p.xlim((-0.7,0.6))
p.ylim((0.5,1.5))
#p.yscale('log')
p.grid()
p.savefig(join(dir,"fit-"+cos+"-differential-function-universal-residual-log-cz.png"))
p.clf()