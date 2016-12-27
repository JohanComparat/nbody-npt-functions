from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import lib_functions_1pt as lib
import cPickle 
import sys
import astropy.cosmology as co
cosmo = co.Planck13

import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

from scipy.optimize import minimize
from scipy.optimize import curve_fit

from scipy.interpolate import interp1d
from scipy.misc import derivative

import scipy.stats as st
#Quantity studied
#=================
qty = "vmax"
# working directory
#=================
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)

# fitting function parameters
#=================
NminCount = 100
limits_04 =  [125, 450] #max : 2e13#[100, 1000]
limits_10 =  [250, 800] #max : 2e14
limits_25 =  [600, 1100] #max : 5e14
limits_40 =  [800, 1400] #max : 2e14

p0 = n.array([-3, 3., 0.3, 1.])

zmin = 1. #-0.01
zmax = 2.3

cos = "sat"

tolerance = 0.06


#=================
# DATA
#=================
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits"))[1].data
# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.vSelection(data, qty, limits_04, limits_10, limits_25,limits_40) 
# minimum number counts selection
nSel = lib.nSelection(data, NminCount, cos)
# altogether
ok = (zSel) & (mSel) & (nSel)
# selection per box :
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')


# x coordinates definition
#=================
vmax = data[qty]
log_vmax = n.log10(vmax)
#print len(vmax), n.min(vmax), n.max(vmax)


# y coordinates
#=================
norm = (100)**3. /(cosmo.H(data["redshift"]).value)**6.
log_VF = n.log10( norm * vmax**3. * data["dNdVdlnV_"+cos])
log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnV_"+cos+"_c"])
#print n.min(log_VF), n.max(log_VF)

# error on y position
#=================
error = (data["std90_pc_"+cos]**2. + data["dN_counts_"+cos]**(-1.))**(0.5)

x_data = log_vmax[ok]
z_data = data['redshift'][ok]
y_data = log_VF[ok]
y_err = error[ok]

"""
# FITTING THE TREND with REDSHIFT
#=======================
param_z0_file=open(join(dir, "vmax-"+cos+"-diff-function-z0-params.pkl"), 'r')
outCF = cPickle.load(param_z0_file)
param_z0_file.close()
A0, v0, a0, b0 = outCF[0]
print "----------------------------------------------------------"
print A0, v0, a0, b0
print "----------------------------------------------------------"
Az = lambda z, A1 : A0 + A1*z**1.# + z**2 *A2
vz = lambda z, v1 : v0 + v1*z**1.
az = lambda z, a1 : a0 + a1*z**1.
bz = lambda z, b1 : b0 + b1*z**1. # +b3*z**3.

ps = [ -1.5, -0.1, 0., 0.1]

vf = lambda v, z, A1,  v1, a1, b1: n.log10( 10**Az(z, A1) *(1.+ (10**v/10**vz(z,  v1))**(-bz(z, b1))) * n.e**(- (10**v/10**vz(z,  v1))**(az(z,a1) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3])#, ps[4])#, ps[5], ps[6])

chi2fun = lambda ps : n.sum( abs(logFun(x_data, z_data, ps) - y_data) / (y_err) )/(len(y_data) - len(ps))

res = minimize(chi2fun, ps, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 50000000000000})
pOpt = res.x
cov = res.direc
outchi2red = chi2fun(pOpt)
outndof=(len(y_data) - len(ps))
print "----------------------------------------------------------"
print "------------------- first order -----------------------"
print "chi2 / ndof = ",outchi2red
print n.round(pOpt,3)
print n.round(abs(cov.diagonal())**0.5,3)
print "----------------------------------------------------------"

outfile=open(join(dir, "vmax-"+cos+"-diff-function-z1-params.pkl"), 'w')
cPickle.dump([pOpt, cov], outfile)
outfile.close()

param_z0_file=open(join(dir, "vmax-"+cos+"-diff-function-z1-params.pkl"), 'r')
outCF = cPickle.load(param_z0_file)
param_z0_file.close()
A1, v1, a1, b1 = outCF[0]

Az = lambda z, A2 : A0 + A1*z**1 + A2*z**2 
vz = lambda z, v2 : v0 + v1*z**1. + v2*z**2 
az = lambda z, a2 : a0 + a1*z**1. + a2*z**2 
bz = lambda z, b2 : b0 + b1*z**1. +b2*z**2.

ps = [ -1.5, -0.1, 0., 0.1]

vf = lambda v, z, A2,  v2, a2, b2: n.log10( 10**Az(z, A2) * (1.+(10**v/10**vz(z,  v2))**(-bz(z, b2))) * n.e**(- (10**v/10**vz(z,  v2))**(az(z,a2) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3])#, ps[4])#, ps[5], ps[6])

chi2fun = lambda ps : n.sum( abs(logFun(x_data, z_data, ps) - y_data) / (y_err) )/(len(y_data) - len(ps))

res = minimize(chi2fun, ps, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 50000000000000})
pOpt = res.x
cov = res.direc
outchi2red = chi2fun(pOpt)
outndof=(len(y_data) - len(ps))
print "----------------------------------------------------------"
print "------------ second order --------------------"
print "chi2 / ndof = ",outchi2red
print n.round(pOpt,3)
print n.round(abs(cov.diagonal())**0.5,3)
print "----------------------------------------------------------"

outfile=open(join(dir, "vmax-"+cos+"-diff-function-z2-params.pkl"), 'w')
cPickle.dump([pOpt, cov], outfile)
outfile.close()


param_z0_file=open(join(dir, "vmax-"+cos+"-diff-function-z2-params.pkl"), 'r')
outCF = cPickle.load(param_z0_file)
param_z0_file.close()
A2, v2, a2, b2 = outCF[0]

Az = lambda z, A3 : A0 + A1*z**1 + A2*z**2  + A3*z**3 
vz = lambda z, v3 : v0 + v1*z**1. + v2*z**2  + v3*z**3 
az = lambda z, a3 : a0 + a1*z**1. + a2*z**2  + a3*z**3 
bz = lambda z, b3 : b0 + b1*z**1. +b2*z**2. + b3*z**3

ps = [ -1.5, -0.1, 0., 0.1]

vf = lambda v, z, A2,  v2, a2, b2: n.log10( 10**Az(z, A2) * (1.+(10**v/10**vz(z,  v2))**(-bz(z, b2))) * n.e**(- (10**v/10**vz(z,  v2))**(az(z,a2) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3])#, ps[4])#, ps[5], ps[6])

chi2fun = lambda ps : n.sum( abs(logFun(x_data, z_data, ps) - y_data) / (y_err) )/(len(y_data) - len(ps))

res = minimize(chi2fun, ps, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 50000000000000})
pOpt = res.x
cov = res.direc
outchi2red = chi2fun(pOpt)
outndof=(len(y_data) - len(ps))

print "----------------------------------------------------------"
print "----------------------third order---------------------"

print "chi2 / ndof = ",outchi2red
print n.round(pOpt,3)
print n.round(abs(cov.diagonal())**0.5,3)
print "----------------------------------------------------------"

outfile=open(join(dir, "vmax-"+cos+"-diff-function-z3-params.pkl"), 'w')
cPickle.dump([pOpt, cov], outfile)
outfile.close()

"""

param_z0_file=open(join(dir, "vmax-"+cos+"-diff-function-z0-params.pkl"), 'r')
outCF = cPickle.load(param_z0_file)
param_z0_file.close()
A0, v0, a0, b0 = outCF[0]
print "----------------------------------------------------------"
print A0, v0, a0, b0
print "----------------------------------------------------------"
"""
Az = lambda z, A1, A3 : A0 + A1*z**1.+ z**3 *A3
vz = lambda z, v1 : v0 + v1*z**1.
az = lambda z, a1 : a0 + a1*z**1.
bz = lambda z, b1, b3 : b0 + b1*z**1. + b3*z**3.

ps = [ -0.6, 0., -0.1, -0.1, 0.4, 0.0]

vf = lambda v, z, A1, A3, v1, a1, b1, b3: n.log10( 10**Az(z, A1, A3) * (1.+ (10**v/10**vz(z,  v1))**(-bz(z, b1, b3))) * n.e**(- (10**v/10**vz(z,  v1))**(az(z, a1) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3], ps[4], ps[5])#, ps[6])

Az = lambda z, A1, A2 : A0 + A1*z**1.+ A2*z**3.
vz = lambda z, v1, v2 : v0 + v1*z**1. + v2*z**3.
az = lambda z, a1, a2 : a0 + a1*z**1. + a2*z**3.
bz = lambda z, b1, b2 : b0 + b1*z**1. + b2*z**3.

ps = [ -0.6, 0., -0.1, 0., -0.1, 0., 0.4, 0.0]

vf = lambda v, z, A1, A2, v1, v2, a1, a2, b1, b2: n.log10( 10**Az(z, A1, A2) * (1.+ (10**v/10**vz(z,  v1, v2))**(-bz(z, b1, b2))) * n.e**(- (10**v/10**vz(z,  v1, v2))**(az(z, a1, a2) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3], ps[4], ps[5], ps[6], ps[7])
"""

Az = lambda z, A0, A1 : A0 + A1*z**1.
vz = lambda z, v0, v1 : v0 + v1*z**1. 
az = lambda z, a0 : a0# + a1*z**1. 
bz = lambda z, b0 : b0# + b1*z**1. 

ps = [ A0, -0.58, v0, -0.2, a0 , b0]

vf = lambda v, z, A0, A1, v0, v1, a0, b0: n.log10( 10**Az(z, A0, A1) * (1.+ (10**v/10**vz(z,  v0, v1))**(-bz(z, b0))) * n.e**(- (10**v/10**vz(z,  v0, v1))**(az(z, a0) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3], ps[4], ps[5])

chi2fun = lambda ps : n.sum( abs(logFun(x_data, z_data, ps) - y_data) / (y_err) )/(len(y_data) - len(ps))
chi2funpp = lambda ps :  abs(logFun(x_data, z_data, ps) - y_data) / (y_err)/(len(y_data) - len(ps))

res = minimize(chi2fun, ps, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 50000000000000})
pOpt = res.x
cov = res.direc
outchi2red = chi2fun(pOpt)
outndof=(len(y_data) - len(ps))
print "----------------------------------------------------------"
print "chi2 / ndof = ", outchi2red, outndof, outndof*outchi2red
print n.round(pOpt,3)
print n.round(abs(cov.diagonal())**0.5,3)
print "----------------------------------------------------------"

outfile=open(join(dir, "vmax-"+cos+"-diff-function-z1-params.pkl"), 'w')
cPickle.dump([pOpt, cov], outfile)
outfile.close()

outfile=open(join(dir, "vmax-"+cos+"-diff-function-z1-params.pkl"), 'w')
cPickle.dump(res, outfile)
outfile.close()

X = n.arange(n.min(x_data),n.max(x_data), 0.01)
Z = n.arange(zmin, zmax, 0.01)
x_model, z_model = n.meshgrid(X, Z)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(x_data, logFun(x_data, z_data, pOpt) , c=z_data, s=5, marker='o',label="model", rasterized=True, vmin=zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'log$_{10} [V^3/H^3(z)\; dn(V)/dlnV]$') # log$_{10}[ n(>M)]')
p.ylim((-5.5,0))
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
#p.yscale('log')
p.grid()
p.savefig(join(dir, "vmax-"+cos+"-differential-function-redshift-trend-model.png"))
p.clf()

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(x_data, chi2funpp(pOpt) , c=z_data, s=5, marker='o',label="model", rasterized=True, vmin=zmin, vmax = zmax)
oneLine = 1./(float(len(x_data))-len(pOpt))
p.axhline(oneLine, label='1/(ndof)')
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'$\chi^2/$point') # log$_{10}[ n(>M)]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.ylim((0.01*oneLine, 100*oneLine))
p.yscale('log')
p.grid()
p.savefig(join(dir, "vmax-"+cos+"-differential-function-redshift-trend-chi2pp.png"))
p.clf()
"""

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(-x_data, y_data , c=z_data, s=5, marker='o',label="data", rasterized=True, vmin=zmin, vmax=zmax)
#p.errorbar(-x_data, y_data , yerr=y_err,fmt='none', rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'log$_{10} V^4 n(>V)$') # log$_{10}[ n(>M)]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.xlim((-0.6, 0.4))
p.ylim((-3, 0.))
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"vmax-"+cos+"-differential-function-redshift-trend-data.png"))
p.clf()

"""

MD_sel_fun=lambda name : (data["boxName"]==name)
MDnames= n.array(['MD_0.4Gpc', 'MD_1Gpc', 'MD_2.5Gpc','MD_4Gpc','MD_2.5GpcNW','MD_4GpcNW'])
MDsels=n.array([MD_sel_fun(name)[ok] for name in MDnames])

f_diff_fun = lambda MDs:  y_data[MDs] - logFun(x_data[MDs], z_data[MDs], pOpt)
f_diffs = n.array([f_diff_fun(MD) for MD in MDsels])

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
for index, fd in enumerate(f_diffs):
	inTol = (abs(10**fd-1)<tolerance)
	print index
	if len(fd)>0:
		p.errorbar(x_data[MDsels[index]], 10**fd, yerr = y_err[MDsels[index]] , rasterized=True, fmt='none', label=MDnames[index])
		print len(inTol.nonzero()[0]), len(fd), 100.*len(inTol.nonzero()[0])/ len(fd)

p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)

p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'data / model') 
#p.xlim((-0.6, 0.4))
#p.yscale('log')
p.grid()
p.savefig(join(dir,"vmax-"+cos+"-differential-function-redshift-trend-residual.png"))
p.clf()

"""
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(x_data, f_diff , c=z_data, s=5, marker='o',label="model reiduals", rasterized=True, vmin=zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)

p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'data - model') 
#p.xlim((-0.6, 0.4))
#p.ylim((-0.9, 1.1))
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"vmax-"+cos+"-differential-function-redshift-trend-residual-z.png"))
p.clf()


tolerance = 0.06
in04 = (abs(10**f_diff_04-1)<tolerance)
print len(in04.nonzero()[0]), len(f_diff_04), 100.*len(in04.nonzero()[0])/ len(f_diff_04)
in10 = (abs(10**f_diff_10-1)<tolerance)
print len(in10.nonzero()[0]), len(f_diff_10), 100.*len(in10.nonzero()[0])/ len(f_diff_10)
in25 = (abs(10**f_diff_25-1)<tolerance)
print len(in25.nonzero()[0]), len(f_diff_25), 100.*len(in25.nonzero()[0])/ len(f_diff_25)
in40 = (abs(10**f_diff_40-1)<tolerance)
print len(in40.nonzero()[0]), len(f_diff_40), 100.*len(in40.nonzero()[0])/ len(f_diff_40)

tolerance = 0.08
inall = (abs(10**f_diff-1)<tolerance)
print len(inall.nonzero()[0]), len(f_diff), 100.*len(inall.nonzero()[0])/ len(f_diff)
"""