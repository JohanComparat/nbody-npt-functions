from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import sys
from scipy import linalg

import lib_functions_1pt as lib

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu

import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p


# mass function theory
from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

#lib.covariance_factor
#lib.f_BH(sigma, 0.333, 0.788, 0.807, 1.795)
A0=0.290
a0=0.8915
p0=0.5524
q0=1.578
#bias = lambda sigma, a0, p0, q0 : lib.b_BH(sigma, a0, p0, q0)
bias = lambda sigma : lib.b_BH(sigma, a=0.8915, p=0.5524, q=1.578)
fsigma = lambda sigma : lib.f_BH(sigma, A0, a0, p0, q0)

#Quantity studied
qty = "mvir"
# working directory
dir = join(os.environ['MVIR_DIR'])
# loads summary file
data = fits.open( join(dir, qty+"_summary.fits"))[1].data

NminCount = 10000
logNpmin = 4.1

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

# We use the 2.5 Gpc/h setup for the covariance matrix to be rescaled to the need
binW=0.025
fc = lib.covariance_factor_jk2[2]
growth_factor_square = 1 # (lib.hmf.growth_factor)**2.
lbox = 400.
snFactor = n.polyval([2.95, -0.04], n.log10(lbox))
svFactor = 2.

var_total_sv = lambda s1, s2 : bias(s1) * bias(s2) * fc / svFactor # * growth_factor_square 
var_total_sn = lambda s1, s2 : lib.shot_double(s1, s2, ( lbox)**3. , binW) / snFactor 
var_total = lambda s1, s2 : (var_total_sn(s1, s2) + var_total_sv(s1, s2))*1e5


sigma_data = 10**-logsig[MD04]
f_data = 10**log_MF[MD04]
f_data_err = (data["std90_pc_cen"][MD04]**2. + data["dN_counts_cen"][MD04]**(-1.))**(0.5)
print len(f_data)

s1, s2 = n.meshgrid(sigma_data, sigma_data)
C = var_total(s1, s2)
#C=n.array([[1.,2.],[3.,5.]])
print n.linalg.cond(C), n.finfo(C.dtype).eps

U, s, Vh = linalg.svd(C, full_matrices=False)
S=n.diag(s)
S_inv=n.diag(s**(-1.))
#print n.dot(S, S_inv)
#print n.allclose(C, n.dot(U, n.dot(S, Vh)))

C_inv_a = n.dot(Vh.T, n.dot(S_inv, U.T))
C_inv_b = n.dot(n.dot(Vh.T, S_inv), U.T)
C_inv_c = n.linalg.inv(C)
print n.dot(C, C_inv_a).diagonal()
print n.dot(C, C_inv_b).diagonal()
print n.dot(C, C_inv_c).diagonal()

C_inv = C_inv_cqey
sys.exit()
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
