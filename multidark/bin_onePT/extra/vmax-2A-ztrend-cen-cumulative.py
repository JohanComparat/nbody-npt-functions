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

data = fits.open( join("..", "vmax", "MD_vmax_summary.fits") )[1].data

errorFactor = 3.
systError = 0.1
NminCount = 1000
limits_04 = [100, 400]
limits_10 = [250, 1000]
limits_25 = [600, 1300]
limits_40 = [1200, 1600]
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin = -0.01
zmax = 2.3
qty = 'vmax'
cos = "cen"

outfile=open(join(dir,qty,"vmax-"+cos+"-cumul-function-z0-params.pkl"), 'r')
outCF = cPickle.load(outfile)
outfile.close()
A0, v0, a0, b0 = outCF[0]
print "----------------------------------------------------------"
print A0, v0, a0, b0
print "----------------------------------------------------------"

"""
print " test 6 params" # 2.79
Az = lambda z, A1, A2 : A0 + z*A1 + A2* z**2
vz = lambda z, v1, v2: v0 + z*v1 + v2*z**2.
az = lambda z : a0 #+ a1*z + a2*z**2
bz = lambda z, b1, b2 : b0 +z*b1 +b2*z**2.
ps = [ 0., 0., 0., 0., 0., 0.]
vf = lambda v, z, A1, A2, v1, v2, b1, b2: n.log10( 10**Az(z, A1, A2) * (10**v/10**vz(z,  v1, v2))**(-bz(z, b1, b2)) * n.e**(- (10**v/10**vz(z, v1, v2))**(az(z) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3], ps[4], ps[5])

print " test 8 params"
Az = lambda z, A1, A2 : A0 + z*A1 + A2* z**2
vz = lambda z, v1, v2: v0 + z*v1 + v2*z**2.
az = lambda z, a1, a2 : a0 + a1*z + a2*z**2
bz = lambda z, b1, b2 : b0 +z*b1 +b2*z**2.
ps = [ 0., 0., 0., 0., 0., 0., 0., 0.]
vf = lambda v, z, A1, A2, v1, v2, a1, a2, b1, b2: n.log10( 10**Az(z, A1, A2) * (10**v/10**vz(z,  v1, v2))**(-bz(z, b1, b2)) * n.e**(- (10**v/10**vz(z, v1, v2))**(az(z, a1, a2) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3], ps[4], ps[5], ps[6], ps[7])

print " test 4 params" # 2.90 -1.45 -0.83 0.096 0.368
Az = lambda z, A1 : A0 + z*A1# + A2* z**2
vz = lambda z, v1: v0 + z*v1# + v2*z**2.
az = lambda z, a1 : a0 + a1*z# + a2*z**2
bz = lambda z, b1 : b0 +z*b1# +b2*z**2.
ps = [ 0., 0., 0., 0.]
vf = lambda v, z, A1, v1, a1, b1: n.log10( 10**Az(z, A1) * (10**v/10**vz(z,  v1))**(-bz(z, b1)) * n.e**(- (10**v/10**vz(z, v1))**(az(z, a1) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3])


print " test 3 params" # 2.96, -1.392, -0.103, 0.314
Az = lambda z, A1 : A0 + z*A1# + A2* z**2
vz = lambda z, v1: v0 + z*v1# + v2*z**2.
az = lambda z : a0 
bz = lambda z, b1 : b0 +z*b1# +b2*z**2.
ps = [ 0., 0., 0.]
vf = lambda v, z, A1, v1, b1: n.log10( 10**Az(z, A1) * (10**v/10**vz(z,  v1))**(-bz(z, b1)) * n.e**(- (10**v/10**vz(z, v1))**(az(z) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2])


print " test 5 params" # 2.78
Az = lambda z, A1, A2 : A0 + z*A1 + A2* z**2
vz = lambda z, v1: v0 + z*v1 #+ v2*z**2.
az = lambda z : a0 #+ a1*z + a2*z**2
bz = lambda z, b1, b2 : b0 +z*b1 +b2*z**3.
ps = [ -1.4, 0., -0.1, 0.27, 0.0]
vf = lambda v, z, A1, A2, v1, b1, b2: n.log10( 10**Az(z, A1, A2) * (10**v/10**vz(z,  v1))**(-bz(z, b1, b2)) * n.e**(- (10**v/10**vz(z, v1))**(az(z) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3], ps[4])
"""
print " test 5 params" # 2.78
Az = lambda z, A1, A2 : A0 + z*A1 + A2* z**2
vz = lambda z, v1: v0 + z*v1 #+ v2*z**2.
az = lambda z : a0 #+ a1*z + a2*z**2
bz = lambda z, b1 : b0 +z*b1# +b2*z**3.
ps = [ -1.4, 0., -0.1, 0.27]
vf = lambda v, z, A1, A2, v1, b1: n.log10( 10**Az(z, A1, A2) * (10**v/10**vz(z,  v1))**(-bz(z, b1)) * n.e**(- (10**v/10**vz(z, v1))**(az(z) ) ) )
logFun = lambda v, z, ps : vf(v, z, ps[0], ps[1], ps[2], ps[3])

# redshift selection
zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
# mass selection
mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>limits_04[0]) &(data["log_"+qty+"_max"]<limits_04[1])) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>limits_10[0]) &(data["log_"+qty+"_max"]<limits_10[1])) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>limits_25[0]) &(data["log_"+qty+"_max"]<limits_25[1])) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>limits_40[0])&(data["log_"+qty+"_max"]<limits_40[1])) 
# minimum number counts selection
nSel = (data['dN_counts_'+cos]>NminCount)
# altogether
ok = (zSel) & (mSel) & (nSel)

MD04=(ok)&(data["boxLength"]==400.)
MD10=(ok)&(data["boxLength"]==1000.)
MD25=(ok)&(data["boxLength"]==2500.)
MD40=(ok)&(data["boxLength"]==4000.)

x_data_04 = (n.log10(data["log_"+qty+"_min"][MD04])+n.log10(data["log_"+qty+"_max"][MD04]))/2.
y_err_04 = (data['std90_pc_cen'][MD04]*errorFactor+systError)/ n.log(10.) 
norm = (100)**3. /(cosmo.H(data["redshift"][MD04]).value)**6.
y_data_04 = n.log10( norm * (10**x_data_04)**3. * data["dNdVdlnM_"+cos+"_c"][MD04])
z_data_04 = data["redshift"][MD04]

x_data_10 = (n.log10(data["log_"+qty+"_min"][MD10])+n.log10(data["log_"+qty+"_max"][MD10]))/2.
y_err_10 = (data['std90_pc_cen'][MD10]*errorFactor+systError)/ n.log(10.)
norm = (100)**3. /(cosmo.H(data["redshift"][MD10]).value)**6.
y_data_10 = n.log10( norm * (10**x_data_10)**3. * data["dNdVdlnM_"+cos+"_c"][MD10]) 
z_data_10 = data["redshift"][MD10]

x_data_25 = (n.log10(data["log_"+qty+"_min"][MD25]) +n.log10(data["log_"+qty+"_max"][MD25]))/2.
y_err_25 = (data['std90_pc_cen'][MD25]*errorFactor+systError)/ n.log(10.) 
norm = (100)**3. /(cosmo.H(data["redshift"][MD25]).value)**6.
y_data_25 = n.log10( norm * (10**x_data_25)**3. * data["dNdVdlnM_"+cos+"_c"][MD25])
z_data_25 = data["redshift"][MD25]

x_data_40 = (n.log10(data["log_"+qty+"_min"][MD40])+n.log10(data["log_"+qty+"_max"][MD40]))/2.
y_err_40 = (data['std90_pc_cen'][MD40]*errorFactor+systError)/ n.log(10.) 
norm = (100)**3. /(cosmo.H(data["redshift"][MD40]).value)**6.
y_data_40 = n.log10( norm * (10**x_data_40)**3. * data["dNdVdlnM_"+cos+"_c"][MD40])
z_data_40 = data["redshift"][MD40]


# x coordinates
x_data = (n.log10(data["log_"+qty+"_min"][ok])+n.log10(data["log_"+qty+"_max"][ok]))/2.
vmax = 10**x_data
# y coordinates
norm = (100)**3. /(cosmo.H(data["redshift"][ok]).value)**6.
y_data = n.log10( norm * (10**x_data)**3. * data["dNdVdlnM_"+cos+"_c"][ok])
y_err = (errorFactor*data['std90_pc_cen'][ok] + systError )/ n.log(10.) 
z_data = data["redshift"][ok]
# chi2fun = lambda ps : n.sum( (logFun(x_data, ps) - y_data)**2. / (y_err)**2. )/(len(y_data) - len(ps))
chi2fun = lambda ps : n.sum( abs(logFun(x_data, z_data, ps) - y_data) / (y_err) )/(len(y_data) - len(ps))



res = minimize(chi2fun, ps, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 50000000000000})
pOpt = res.x
cov = res.direc
#chi2perpoint = lambda ps : (funG(lg_vmax, lg_1pz, ps) - lg_MF_c)**2. / (errorLog)**2. 
#chi2pp = chi2perpoint(pOpt)
#|print pOpt, cov
print "----------------------------------------------------------"
print n.round(pOpt,4)
print n.round(abs(cov.diagonal())**0.5,4)
print "----------------------------------------------------------"

outfile=open(join(dir,qty,"vmax-"+cos+"-cumul-function-z1-params.pkl"), 'w')
cPickle.dump(res, outfile)
outfile.close()

X = n.arange(n.min(x_data),n.max(x_data), 0.01)
Z = n.arange(zmin, zmax, 0.01)
x_model, z_model = n.meshgrid(X, Z)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])

#sc1=p.scatter(-n.hstack((x_model)), logFun(n.hstack((x_model)), n.hstack((z_model)), pOpt) , c=n.hstack((z_model)), s=2, marker='o',label="model", rasterized=True, vmin=zmin, vmax=zmax)

sc1=p.scatter(x_data, logFun(x_data, z_data, pOpt) , c=z_data, s=5, marker='o',label="model", rasterized=True, vmin=zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'log$_{10}( V^3/H^3 n(>V) )$') # log$_{10}[ n(>M)]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"vmax-"+cos+"-cumulative-function-redshift-trend-model.png"))
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
p.savefig(join(dir,qty,"vmax-"+cos+"-cumulative-function-redshift-trend-data.png"))
p.clf()

"""
f_diff_04 = y_data_04 - logFun(x_data_04, z_data_04, pOpt) 
f_diff_10 = y_data_10 - logFun(x_data_10, z_data_10, pOpt) 
f_diff_25 = y_data_25 - logFun(x_data_25, z_data_25, pOpt) 
f_diff_40 = y_data_40 - logFun(x_data_40, z_data_40, pOpt) 
f_diff = y_data - logFun(x_data, z_data, pOpt) 

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.errorbar(x_data_04, 10**f_diff_04, yerr = y_err_04 , rasterized=True, fmt='none', label="MD04")
p.errorbar(x_data_10, 10**f_diff_10, yerr = y_err_10 , rasterized=True, fmt='none', label="MD10")
p.errorbar(x_data_25, 10**f_diff_25, yerr = y_err_25 , rasterized=True, fmt='none', label="MD25")
p.errorbar(x_data_40, 10**f_diff_40, yerr = y_err_40 , rasterized=True, fmt='none', label="MD40")
p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)

p.xlabel(r'log$_{10}[V_{max}/(km \; s^{-1})]$')
p.ylabel(r'data / model') 
#p.xlim((-0.6, 0.4))
#p.ylim((-0.9, 1.1))
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"vmax-"+cos+"-cumulative-function-redshift-trend-residual.png"))
p.clf()


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
p.savefig(join(dir,qty,"vmax-"+cos+"-cumulative-function-redshift-trend-residual-z.png"))
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
