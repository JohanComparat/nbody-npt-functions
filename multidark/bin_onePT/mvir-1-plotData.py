from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import sys

import lib_functions_1pt as lib

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu

import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p

#Quantity studied
qty = "mvir"
# working directory
dir = join(os.environ['MVIR_DIR'])
# loads summary file
data = fits.open( join(dir, qty+"_summary.fits"))[1].data

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

# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, logNpmin)
mSel2_inter = (data["log_mvir"]<13.2) & (data["redshift"]>0.)
mSel2 = (mSel2_inter==False)
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (mSel2) & (nSelCen)
# selection per box :
MD04=(ok)&(data["boxName"]=='MD_0.4Gpc')
MD10=(ok)&(data["boxName"]=='MD_1Gpc_new_rockS')
MD25=(ok)&(data["boxName"]=='MD_2.5Gpc')
MD40=(ok)&(data["boxName"]=='MD_4Gpc')
MD25NW=(ok)&(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(ok)&(data["boxName"]=='MD_4GpcNW')
DS80=(ok)&(data["boxName"]=='DS_8Gpc')

ok = (zSel) & (mSel) & (mSel2) & (nSelCen)&(data["boxName"]!='DS_8Gpc')

# NOW PLOTTING ALL THE DATA
#lib.plot_mvir_function_data(log_mvir[ok], logsig[ok], lognu[ok], log_MF[ok], log_MF_c[ok], data['redshift'][ok], zmin, zmax, cos = cos)

# ERROR PLOT: JK vs. POISSON
x = data["std90_pc_"+cos] 
y = data["dN_counts_"+cos]**(-0.5)
#lib.plot_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, DS80, cos = cos, dir=join(os.environ['MVIR_DIR']))




p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
# log_mvir, logsigM1, logNu, log_MF, log_MF_c, redshift  = log_mvir[ok], logsig[ok], lognu[ok], log_MF[ok], log_MF_c[ok], data['redshift'][ok]
x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_cen"][ok]**2. + data["dN_counts_cen"][ok]**(-1.))**(0.5)
p.errorbar(x_data, y_data, yerr = y_data_err, rasterized=True, fmt='none', label='distinct halos z=0', lw=2)
p.plot(lib.X, n.log10(lib.ftC16), 'k--', label='fit', lw=2)

cos = 'sat'
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )
nSelSat = lib.nSelection(data, NminCount, cos )
ok = (zSel) & (mSel)& (mSel2)  & (nSelSat)&(data["boxName"]!='DS_8Gpc')

x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_"+cos][ok]**2. + data["dN_counts_"+cos][ok]**(-1.))**(0.5)
p.errorbar(x_data, y_data, yerr = y_data_err, rasterized=True, fmt='none', label='satellite subhalos z=0', lw=2)
sigs = n.arange(-0.5,.6, 0.01)
p.plot(lib.X, n.log10(lib.ftC16st_sat), 'k--', lw=2)

p.xlabel(r'$log_{10}(\sigma^{-1})$')
p.ylabel(r'$\log_{10}\left[ \frac{M}{\rho_m} \frac{dn}{d\ln M} \left|\frac{d\ln M }{d\ln \sigma}\right|\right] $') 
 # log$_{10}[ n(>M)]')
gl = p.legend(loc=0, fontsize=12)
gl.set_frame_on(False)
p.ylim((-3., -0.4))
p.xlim((-0.5, 0.5))
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],"mvir-AL-z0-differential-function-data-xSigma.png"))
p.clf()

sys.exit()
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

# minimum number counts selection
nSelSat = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel)& (mSel2)  & (nSelSat)

# NOW PLOTTING ALL THE DATA
lib.plot_mvir_function_data(log_mvir[ok], logsig[ok], lognu[ok], log_MF[ok], log_MF_c[ok], data['redshift'][ok], zmin, zmax, cos = cos)

lib.plot_mvir_function_data_perBox(log_mvir, log_MF, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos=cos)

# ERROR PLOT: JK vs. POISSON
x = data["std90_pc_"+cos] 
y = data["dN_counts_"+cos]**(-0.5)
lib.plot_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, DS80, cos = cos, dir=join(os.environ['MVIR_DIR']))

