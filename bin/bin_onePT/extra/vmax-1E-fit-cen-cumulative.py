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

# fitting function
p0 = n.array([-1.7, 3., 1.8, 0.8])
vf = lambda v, A, v0, alpha, beta : n.log10( 10**A * (10**v/10**v0)**(-beta) * n.e**(- (10**v/10**v0)**(alpha) ) )

errorFactor = 3.
systError = 0.01

dir='..'
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")
dir_25N = join(dir,"MD_2.5GpcNW")
dir_40N = join(dir,"MD_4GpcNW")

data = fits.open( join("..", "vmax", "MD_vmax_summary.fits") )[1].data

NminCount = 1000
limits_04 = [100, 400]
limits_10 = [250, 1000]
limits_25 = [600, 1300]
limits_40 = [1200, 1600]
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin =  -0.01
zmax = 0.001
qty = 'vmax'
cos = "cen"

# redshift selection
zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
# mass selection
if  cos == "cen":
	mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>limits_04[0]) &(data["log_"+qty+"_max"]<limits_04[1])) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>limits_10[0]) &(data["log_"+qty+"_max"]<limits_10[1])) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>limits_25[0]) &(data["log_"+qty+"_max"]<limits_25[1])) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>limits_40[0])&(data["log_"+qty+"_max"]<limits_40[1])) 
if  cos == "sat":
	mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>limits_04[0]) &(data["log_"+qty+"_max"]<limits_04[1])) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>limits_10[0]) &(data["log_"+qty+"_max"]<limits_10[1])) 

# minimum number counts selection
nSel = (data['dN_counts_'+cos]>NminCount)
# altogether
ok = (zSel) & (mSel) & (nSel)
	
MD04=(ok)&(data["boxLength"]==400.)
MD10=(ok)&(data["boxLength"]==1000.)
MD25=(ok)&(data["boxLength"]==2500.)
MD40=(ok)&(data["boxLength"]==4000.)

x_data_04 = (n.log10(data["log_"+qty+"_min"][MD04])+n.log10(data["log_"+qty+"_max"][MD04]))/2.
error_04 = (data['std90_pc_cen'][MD04]*errorFactor+systError)/ n.log(10.) 
norm = (100)**3. /(cosmo.H(data["redshift"][MD04]).value)**6.
y_data_04 = n.log10( norm * (10**x_data_04)**3. * data["dNdVdlnM_"+cos+"_c"][MD04])
#y_data_04 = n.log10( (10**x_data_04)**4. * data["dNdVdlnM_"+cos][MD04])

x_data_10 = (n.log10(data["log_"+qty+"_min"][MD10])+n.log10(data["log_"+qty+"_max"][MD10]))/2.
error_10 = (data['std90_pc_cen'][MD10]*errorFactor+systError)/ n.log(10.)
norm = (100)**3. /(cosmo.H(data["redshift"][MD10]).value)**6.
y_data_10 = n.log10( norm * (10**x_data_10)**3. * data["dNdVdlnM_"+cos+"_c"][MD10]) 
#y_data_10 = n.log10( (10**x_data_10)**4. * data["dNdVdlnM_"+cos][MD10])

x_data_25 = (n.log10(data["log_"+qty+"_min"][MD25]) +n.log10(data["log_"+qty+"_max"][MD25]))/2.
error_25 = (data['std90_pc_cen'][MD25]*errorFactor+systError)/ n.log(10.) 
norm = (100)**3. /(cosmo.H(data["redshift"][MD25]).value)**6.
y_data_25 = n.log10( norm * (10**x_data_25)**3. * data["dNdVdlnM_"+cos+"_c"][MD25])
#y_data_25 = n.log10( (10**x_data_25)**4. * data["dNdVdlnM_"+cos][MD25])

x_data_40 = (n.log10(data["log_"+qty+"_min"][MD40])+n.log10(data["log_"+qty+"_max"][MD40]))/2.
error_40 = (data['std90_pc_cen'][MD40]*errorFactor+systError)/ n.log(10.) 
norm = (100)**3. /(cosmo.H(data["redshift"][MD40]).value)**6.
y_data_40 = n.log10( norm * (10**x_data_40)**3. * data["dNdVdlnM_"+cos+"_c"][MD40])
#y_data_40 = n.log10( (10**x_data_40)**4. * data["dNdVdlnM_"+cos][MD40])

# x coordinates
x_data = (n.log10(data["log_"+qty+"_min"][ok])+n.log10(data["log_"+qty+"_max"][ok]))/2.
vmax = 10**x_data
# y coordinates
norm = (100)**3. /(cosmo.H(data["redshift"][ok]).value)**6.
y_data = n.log10( norm * (10**x_data)**3. * data["dNdVdlnM_"+cos+"_c"][ok])
y_err = (errorFactor*data['std90_pc_cen'][ok] + systError )/ n.log(10.) 
	
outCF=curve_fit(vf, x_data, y_data, p0, y_err)#, bounds=boundaries)

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
#chi2perpoint = lambda ps : (funG(lg_vmax, lg_1pz, ps) - lg_MF_c)**2. / (errorLog)**2. 
#chi2pp = chi2perpoint(pOpt)
print pOpt, cov
"""
x_model = n.arange(n.min(x_data),n.max(x_data),0.005)
y_model = vf(x_model, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
n.savetxt(join(dir,qty,"vmax-"+cos+"-cumulative-function-z0-model-pts.txt"),n.transpose([x_model, y_model]) )

outfile=open(join(dir,qty,"vmax-"+cos+"-cumul-function-z0-params.pkl"), 'w')
cPickle.dump(outCF, outfile)
outfile.close()



f_diff =  y_data - vf(x_data, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
print "chi2 red=", n.sum((f_diff/y_err)**2.)/len(y_err-4), len(y_err-4)
print "params=",outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3]
print "err=",outCF[1][0][0]**0.5, outCF[1][1][1]**0.5, outCF[1][2][2]**0.5, outCF[1][3][3]**0.5
f_diff_04 =  y_data_04 - vf(x_data_04, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
f_diff_40 =  y_data_40 - vf(x_data_40, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
f_diff_10 =  y_data_10 - vf(x_data_10, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])
f_diff_25 =  y_data_25 - vf(x_data_25, outCF[0][0], outCF[0][1], outCF[0][2], outCF[0][3])


# now the plots
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.errorbar(x_data_04, 10**f_diff_04, yerr = error_04 , rasterized=True, fmt='none', label="MD04")
p.errorbar(x_data_10, 10**f_diff_10, yerr = error_10 , rasterized=True, fmt='none', label="MD10")
p.errorbar(x_data_25, 10**f_diff_25, yerr = error_25 , rasterized=True, fmt='none', label="MD25")
p.errorbar(x_data_40, 10**f_diff_40, yerr = error_40 , rasterized=True, fmt='none', label="MD40")
p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
p.axhline(0.99,c='k',ls='--')
p.xlabel(r'$log(V_{max})$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
#p.xlim((-0.7,0.6))
#p.ylim((-0.05,0.05))
#p.yscale('log')
p.grid()
#p.title(n.round(outCF[0],3))
p.savefig(join(dir,qty,"vmax-"+cos+"-cumulative-function-fit-residual-log.png"))
p.clf()


tolerance = 0.03
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

