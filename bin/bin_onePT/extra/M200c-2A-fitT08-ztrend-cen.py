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

dir='..'
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")
dir_25N = join(dir,"MD_2.5GpcNW")
dir_40N = join(dir,"MD_4GpcNW")

data = fits.open( join("..", "M200c", "MD_M200c_summary.fits") )[1].data

errorFactor = 3.
systError = 0.01
NminCount = 1000
Npmin = 1000
limits_04 = [Npmin*9.63 * 10**7, 5e12]
limits_10 = [Npmin*1.51 * 10**9., 5e13]
limits_25 = [Npmin*2.359 * 10**10., 5e14]
limits_40 = [Npmin* 9.6 * 10**10. , 5e15]
MPART = n.array([9.63 * 10**7, 1.51 * 10**9, 2.359 * 10**10, 9.6 * 10**10])
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin = -0.01
zmax = 2.3
qty = 'M200c'
cos = "cen"

outfile=open(join(dir,qty,"M200c-"+cos+"-diff-function-z0-params.pkl"), 'r')
outCF = cPickle.load(outfile)
outfile.close()
A0, a0, b0, c0 = outCF[0]
Az = lambda z, A1 : A0*(1+z)**A1
az = lambda z, a1 : a0*(1+z)**a1
bz = lambda z, b1 : b0*(1+z)**b1
#cz = lambda z, c1 : c0#*(1+z)**c1

ps = [-0.1, -0.1, 0.]#, -0.2]

sigma = n.arange(0.05,10,0.05)
f = lambda sigma, z, A1, a1, b1 : Az(z, A1)*( (sigma/bz(z, b1))**(-az(z, a1)) + 1 )*n.e**(-c0/sigma**2.)

logFun = lambda logSigma, z, ps : n.log10( f(10.**logSigma, z, ps[0], ps[1], ps[2]) )

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
#print n.min(logsigM1), n.max(logsigM1)
log_m200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
m200c = 10**log_m200c
# mean density array normalization
rhom = cosmo.critical_density(data["redshift"][ok]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][ok])/(100*uu.km/(uu.Mpc*uu.s)))**1.

x_data = n.log10(data['sigmaM'][ok])#
z_data = data['redshift'][ok]#
f_data = (m200c * data["dNdVdlnM_"+cos][ok]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"][ok]) )
y_data = n.log10(f_data)
y_err = (errorFactor*data['std90_pc_cen'][ok] + systError )/ n.log(10.) 

MD04=(ok) & (data["boxLength"]==400.)
MD10=(ok) & (data["boxLength"]==1000.)
MD25=(ok) & (data["boxLength"]==2500.)
MD40=(ok) & (data["boxLength"]==4000.)

z_data_04 = data['redshift'][MD04]#
rhom_04 = cosmo.critical_density(data["redshift"][MD04]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][MD04])/(100*uu.km/(uu.Mpc*uu.s)))**1.
m200c_04 = 10**((data["log_"+qty+"_min"][MD04]+data["log_"+qty+"_max"][MD04])/2.)
x_data_04 = n.log10(data['sigmaM'][MD04])
y_data_04 = n.log10((m200c_04 * data["dNdVdlnM_"+cos][MD04]/ rhom_04.value  / abs(data["dlnsigmaM1_o_dlnM"][MD04]) ))
y_err_04 = (errorFactor*data['std90_pc_cen'][MD04] + systError ) / n.log(10.)

z_data_10 = data['redshift'][MD10]#
rhom_10 = cosmo.critical_density(data["redshift"][MD10]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][MD10])/(100*uu.km/(uu.Mpc*uu.s)))**1.
m200c_10 = 10**((data["log_"+qty+"_min"][MD10]+data["log_"+qty+"_max"][MD10])/2.)
x_data_10 = n.log10(data['sigmaM'][MD10])
y_data_10 = n.log10((m200c_10 * data["dNdVdlnM_"+cos][MD10]/ rhom_10.value  / abs(data["dlnsigmaM1_o_dlnM"][MD10]) ))
y_err_10 =( errorFactor*data['std90_pc_cen'][MD10]  + systError )/ n.log(10.)

z_data_25 = data['redshift'][MD25]#
rhom_25 = cosmo.critical_density(data["redshift"][MD25]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][MD25])/(100*uu.km/(uu.Mpc*uu.s)))**1.
m200c_25 = 10**((data["log_"+qty+"_min"][MD25]+data["log_"+qty+"_max"][MD25])/2.)
x_data_25 = n.log10(data['sigmaM'][MD25])
y_data_25 = n.log10((m200c_25 * data["dNdVdlnM_"+cos][MD25]/ rhom_25.value  / abs(data["dlnsigmaM1_o_dlnM"][MD25]) ))
y_err_25 = (errorFactor*data['std90_pc_cen'][MD25] + systError ) / n.log(10.)

z_data_40 = data['redshift'][MD40]#
rhom_40 = cosmo.critical_density(data["redshift"][MD40]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][MD40])/(100*uu.km/(uu.Mpc*uu.s)))**1.
m200c_40 = 10**((data["log_"+qty+"_min"][MD40]+data["log_"+qty+"_max"][MD40])/2.)
x_data_40 = n.log10(data['sigmaM'][MD40])
y_data_40 = n.log10((m200c_40 * data["dNdVdlnM_"+cos][MD40]/ rhom_40.value  / abs(data["dlnsigmaM1_o_dlnM"][MD40]) ))
y_err_40 =( errorFactor*data['std90_pc_cen'][MD40]  + systError )/ n.log(10.)

# chi2fun = lambda ps : n.sum( (logFun(x_data, ps) - y_data)**2. / (y_err)**2. )/(len(y_data) - len(ps))
chi2fun = lambda ps : n.sum( abs(logFun(x_data, z_data, ps) - y_data) / (y_err) )/(len(y_data) - len(ps))

res = minimize(chi2fun, ps, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt = res.x
cov = res.direc
#chi2perpoint = lambda ps : (funG(lg_M200c, lg_1pz, ps) - lg_MF_c)**2. / (errorLog)**2. 
#chi2pp = chi2perpoint(pOpt)
#|print pOpt, cov
print n.round(pOpt,2)
print n.round(abs(cov.diagonal())**0.5,2)

outfile=open(join(dir,qty,"M200c-"+cos+"-diff-function-z1-params.pkl"), 'w')
cPickle.dump(res, outfile)
outfile.close()

X = n.arange(-0.4, 0.6, 0.01)
Z = n.arange(zmin, zmax, 0.01)
x_model, z_model = n.meshgrid(X, Z)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])

#sc1=p.scatter(-n.hstack((x_model)), logFun(n.hstack((x_model)), n.hstack((z_model)), pOpt) , c=n.hstack((z_model)), s=2, marker='o',label="model", rasterized=True, vmin=zmin, vmax=zmax)

sc1=p.scatter(-x_data, logFun(x_data, z_data, pOpt) , c=z_data, s=5, marker='o',label="model", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'$log(\sigma^{-1})$')
p.ylabel(r'$\log_{10}(f(\sigma))$') 
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.xlim((-0.6, 0.4))
p.ylim((-3, 0.))
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"M200c-"+cos+"-differential-function-redshift-trend.png"))
p.clf()


p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(-x_data, y_data , c=z_data, s=5, marker='o',label="data", rasterized=True, vmin=zmin, vmax=zmax)
#p.errorbar(-x_data, y_data , yerr=y_err,fmt='none', rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'$log(\sigma^{-1})$')
p.ylabel(r'$\log_{10}(f(\sigma))$') 
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.xlim((-0.6, 0.4))
p.ylim((-3, 0.))
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"M200c-"+cos+"-differential-function-redshift-trend-data.png"))
p.clf()


f_diff_04 = y_data_04 - logFun(x_data_04, z_data_04, pOpt) 
f_diff_10 = y_data_10 - logFun(x_data_10, z_data_10, pOpt) 
f_diff_25 = y_data_25 - logFun(x_data_25, z_data_25, pOpt) 
f_diff_40 = y_data_40 - logFun(x_data_40, z_data_40, pOpt) 
f_diff = y_data - logFun(x_data, z_data, pOpt) 

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.errorbar(-x_data_04, 10**f_diff_04, yerr = y_err_04*n.log(10.) , rasterized=True, fmt='none', label="MD04")
p.errorbar(-x_data_10, 10**f_diff_10, yerr = y_err_10*n.log(10.) , rasterized=True, fmt='none', label="MD10")
p.errorbar(-x_data_25, 10**f_diff_25, yerr = y_err_25*n.log(10.) , rasterized=True, fmt='none', label="MD25")
p.errorbar(-x_data_40, 10**f_diff_40, yerr = y_err_40*n.log(10.) , rasterized=True, fmt='none', label="MD40")
p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)

p.xlabel(r'$log(\sigma^{-1})$')
p.ylabel(r'$f(\sigma)$, data / model') 
#p.xlim((-0.6, 0.4))
#p.ylim((-0.9, 1.1))
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"M200c-"+cos+"-differential-function-redshift-trend-residual.png"))
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
