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

sigma = n.arange(0.05,10,0.05)
f = lambda sigma, A, a, b, c : A*( (sigma/b)**(-a) + 1 )*n.e**(-c/sigma**2.)
logFun = lambda logSigma, ps : n.log10( f(10.**logSigma, ps[0], ps[1], ps[2], ps[3]) )

ft08Y = f(sigma, 0.224, 1.67, 1.80, 1.48) 
ft08X = n.log10(1./sigma)

# mass and mass function from Klypin 16
k16x, k16y=n.loadtxt(join("..", "M200c","previousData","Klypin_2016_MF_z0.txt"),  unpack=True)

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
zmax = 0.001
qty = 'M200c'
cos = "cen"
p0 = [0.224, 1.67, 1.80, 1.48]
p0 = [0.1, 1.2, 3.80, 1.3]
boundaries = [(0.01,0.1,0.1,0.1),(0.12,4,5,4)]
ps = p0
#def fitDataTinker08(qty = 'M200c', cos = "cen", zmin = zmin, zmax = zmax, p0 = [0.224, 1.67, 1.80, 1.48]):
"""
Plots the data to be used in the fits later in the analysis.
"""
# redshift selection
zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
# mass selection
if  cos == "cen":
	mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0])) &(data["log_"+qty+"_max"]<n.log10(limits_04[1]))) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])) &(data["log_"+qty+"_max"]<n.log10(limits_10[1]))) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>n.log10(limits_25[0])) &(data["log_"+qty+"_max"]<n.log10(limits_25[1]))) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>n.log10(limits_40[0]))&(data["log_"+qty+"_max"]<n.log10(limits_40[1]))) 
if  cos == "sat":
	mSel =  ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0])) &(data["log_"+qty+"_max"]<n.log10(limits_04[1]))) |  ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])) &(data["log_"+qty+"_max"]<n.log10(limits_10[1])))
	
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
f_data = (m200c * data["dNdVdlnM_"+cos][ok]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"][ok]) )
y_data = n.log10(f_data)
y_err = (errorFactor*data['std90_pc_cen'][ok] + systError )/ n.log(10.) 

MD04=(ok) & (data["boxLength"]==400.)
MD10=(ok) & (data["boxLength"]==1000.)
MD25=(ok) & (data["boxLength"]==2500.)
MD40=(ok) & (data["boxLength"]==4000.)

rhom_04 = cosmo.critical_density(data["redshift"][MD04]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][MD04])/(100*uu.km/(uu.Mpc*uu.s)))**1.
m200c_04 = 10**((data["log_"+qty+"_min"][MD04]+data["log_"+qty+"_max"][MD04])/2.)
x_data_04 = n.log10(data['sigmaM'][MD04])
y_data_04 = n.log10((m200c_04 * data["dNdVdlnM_"+cos][MD04]/ rhom_04.value  / abs(data["dlnsigmaM1_o_dlnM"][MD04]) ))
y_err_04 = (errorFactor*data['std90_pc_cen'][MD04] + systError ) / n.log(10.)

rhom_10 = cosmo.critical_density(data["redshift"][MD10]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][MD10])/(100*uu.km/(uu.Mpc*uu.s)))**1.
m200c_10 = 10**((data["log_"+qty+"_min"][MD10]+data["log_"+qty+"_max"][MD10])/2.)
x_data_10 = n.log10(data['sigmaM'][MD10])
y_data_10 = n.log10((m200c_10 * data["dNdVdlnM_"+cos][MD10]/ rhom_10.value  / abs(data["dlnsigmaM1_o_dlnM"][MD10]) ))
y_err_10 =( errorFactor*data['std90_pc_cen'][MD10]  + systError )/ n.log(10.)

rhom_25 = cosmo.critical_density(data["redshift"][MD25]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][MD25])/(100*uu.km/(uu.Mpc*uu.s)))**1.
m200c_25 = 10**((data["log_"+qty+"_min"][MD25]+data["log_"+qty+"_max"][MD25])/2.)
x_data_25 = n.log10(data['sigmaM'][MD25])
y_data_25 = n.log10((m200c_25 * data["dNdVdlnM_"+cos][MD25]/ rhom_25.value  / abs(data["dlnsigmaM1_o_dlnM"][MD25]) ))
y_err_25 = (errorFactor*data['std90_pc_cen'][MD25] + systError ) / n.log(10.)

rhom_40 = cosmo.critical_density(data["redshift"][MD40]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][MD40])/(100*uu.km/(uu.Mpc*uu.s)))**1.
m200c_40 = 10**((data["log_"+qty+"_min"][MD40]+data["log_"+qty+"_max"][MD40])/2.)
x_data_40 = n.log10(data['sigmaM'][MD40])
y_data_40 = n.log10((m200c_40 * data["dNdVdlnM_"+cos][MD40]/ rhom_40.value  / abs(data["dlnsigmaM1_o_dlnM"][MD40]) ))
y_err_40 =( errorFactor*data['std90_pc_cen'][MD40]  + systError )/ n.log(10.)


f = lambda sigma, A, a, b, c : A*( (sigma/b)**(-a) + 1 )*n.e**(-c/sigma**2.)
logFun = lambda logSigma, ps : n.log10( f(10.**logSigma, ps[0], ps[1], ps[2], ps[3]) )
logf = lambda logSigma, A, a, b, c  : n.log10( f(10.**logSigma, A, a, b, c ) )

outCF=curve_fit(logf, x_data, y_data, p0, y_err, bounds=boundaries)

#p.errorbar(x_data, y_data, yerr = y_err)
#p.plot(x_data, logf(x_data, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3]), 'k+')
#p.plot(x_data, logf(x_data, p0[0], p0[1], p0[2], p0[3]), 'rx')
#p.show()
"""
# chi2fun = lambda ps : n.sum( (logFun(x_data, ps) - y_data)**2. / (y_err)**2. )/(len(y_data) - len(ps))
chi2fun = lambda ps : n.sum( abs(logFun(x_data, ps) - y_data) / (y_err) )/(len(y_data) - len(ps))

res = minimize(chi2fun, p0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt = res.x
cov = res.direc
#chi2perpoint = lambda ps : (funG(lg_M200c, lg_1pz, ps) - lg_MF_c)**2. / (errorLog)**2. 
#chi2pp = chi2perpoint(pOpt)
print pOpt, cov
"""
x_model = n.arange(n.min(x_data),n.max(x_data),0.005)
y_model = logf(x_model, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
n.savetxt(join(dir,qty,"M200c-"+cos+"-differential-function-z0-model-pts.txt"),n.transpose([x_model, y_model]) )

outfile=open(join(dir,qty,"M200c-"+cos+"-diff-function-z0-params.pkl"), 'w')
cPickle.dump(outCF, outfile)
outfile.close()

f_diff =  y_data - logf(x_data, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
print "chi2 red=", n.sum(f_diff/y_err)/len(y_err-4)
print "params=",outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3]
print "err=",outCF[1][0][0]**0.5, outCF[1][1][1]**0.5, outCF[1][2][2]**0.5, outCF[1][3][3]**0.5
f_diff_04 =  y_data_04 - logf(x_data_04, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
f_diff_40 =  y_data_40 - logf(x_data_40, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
f_diff_10 =  y_data_10 - logf(x_data_10, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
f_diff_25 =  y_data_25 - logf(x_data_25, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])

# now the plots
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.errorbar(x_data_04, f_diff_04, yerr = y_err_04 , rasterized=True, fmt='none', label="MD04")
p.errorbar(x_data_10, f_diff_10, yerr = y_err_10 , rasterized=True, fmt='none', label="MD10")
p.errorbar(x_data_25, f_diff_25, yerr = y_err_25 , rasterized=True, fmt='none', label="MD25")
p.errorbar(x_data_40, f_diff_40, yerr = y_err_40 , rasterized=True, fmt='none', label="MD40")
p.xlabel(r'$log(\sigma)$')
p.ylabel(r'log data - log model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
#p.xlim((-0.7,0.6))
#p.ylim((-0.05,0.05))
#p.yscale('log')
p.grid()
#p.title(n.round(outCF[0],3))
p.savefig(join(dir,qty,"M200c-"+cos+"-differential-function-tinkerFit-residual-log.png"))
p.clf()

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.errorbar(x_data_04, f_diff_04, yerr = y_err_04 , rasterized=True, fmt='none', label="MD04")
p.errorbar(x_data_10, f_diff_10, yerr = y_err_10 , rasterized=True, fmt='none', label="MD10")
p.errorbar(x_data_25, f_diff_25, yerr = y_err_25 , rasterized=True, fmt='none', label="MD25")
p.errorbar(x_data_40, f_diff_40, yerr = y_err_40 , rasterized=True, fmt='none', label="MD40")
p.xlabel(r'$log(\sigma)$')
p.ylabel(r'log data - log model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
#p.xlim((-0.7,0.6))
#p.ylim((-0.05,0.05))
#p.yscale('log')
p.grid()
#p.title(n.round(outCF[0],3))
p.savefig(join(dir,qty,"M200c-"+cos+"-differential-function-tinkerFit-residual-log.png"))
p.clf()
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.errorbar(m200c_04, 10**f_diff_04, yerr = y_err_04*n.log(10.) , rasterized=True, fmt='none', label="MD04")
p.errorbar(m200c_10, 10**f_diff_10, yerr = y_err_10*n.log(10.) , rasterized=True, fmt='none', label="MD10")
p.errorbar(m200c_25, 10**f_diff_25, yerr = y_err_25*n.log(10.) , rasterized=True, fmt='none', label="MD25")
p.errorbar(m200c_40, 10**f_diff_40, yerr = y_err_40*n.log(10.) , rasterized=True, fmt='none', label="MD40")
p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
p.xlabel(r'$log(M_{200c})$')
p.ylabel(r'data / model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
#p.xlim((-0.7,0.6))
#p.ylim((-0.05,0.05))
p.xscale('log')
p.grid()
#p.title(n.round(outCF[0],3))
p.savefig(join(dir,qty,"M200c-"+cos+"-differential-function-tinkerFit-residual-xmass.png"))
p.clf()

tolerance = 0.05
in04 = (abs(10**f_diff_04-1)<tolerance)
print len(in04.nonzero()[0]), len(f_diff_04), 100.*len(in04.nonzero()[0])/ len(f_diff_04)
in10 = (abs(10**f_diff_10-1)<tolerance)
print len(in10.nonzero()[0]), len(f_diff_10), 100.*len(in10.nonzero()[0])/ len(f_diff_10)
in25 = (abs(10**f_diff_25-1)<tolerance)
print len(in25.nonzero()[0]), len(f_diff_25), 100.*len(in25.nonzero()[0])/ len(f_diff_25)
in40 = (abs(10**f_diff_40-1)<tolerance)
print len(in40.nonzero()[0]), len(f_diff_40), 100.*len(in40.nonzero()[0])/ len(f_diff_40)
inall = (abs(10**f_diff-1)<tolerance)
print len(inall.nonzero()[0]), len(f_diff), 100.*len(inall.nonzero()[0])/ len(f_diff)
