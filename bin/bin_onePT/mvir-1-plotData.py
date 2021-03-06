from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import sys

#import lib_functions_1pt as lib

delta_c = 1.686

f_BH = lambda sigma, A, a, p, q: A* (2./n.pi)**(0.5) * ( 1 + (sigma**2./(a*delta_c**2.))**(p) )*(delta_c*a**0.5/sigma)**(q)*n.e**(-a*delta_c**2./(2.*sigma**2.))

b_BH = lambda sigma, a, p, q:  1 + (a*(delta_c/sigma)**2. - q) / delta_c + (2*p/delta_c)/(1 + (a*(delta_c/sigma)**2.)**p)
f_ST = lambda sigma, A, a, p: A* (2./n.pi)**(0.5) * ( 1 + (sigma**2./(a**delta_c*2.))**(p) )*(delta_c*a**0.5/sigma)*n.e**(-a*delta_c**2./(2.*sigma**2.))

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu

import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p

from scipy.interpolate import interp1d


from colossus.cosmology import cosmology
cosmology.setCosmology('planck15')
params = {'flat': True, 'H0': 67.77, 'Om0': 0.307115, 'Ob0': 0.048206, 'sigma8': 0.8228, 'ns': 0.9600, 'relspecies': False}
cosmology.setCosmology('myCosmo', params)
from colossus.lss import peaks
delta_c = peaks.collapseOverdensity(corrections = True, z = 0)
#delta_c = 1.686
nu = n.arange(0.2, 8.0, 0.01)
Masses = peaks.massFromPeakHeight(nu, 0.0)
sigma = ( delta_c / nu)



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

tolerance = 0.03


mSelection = lambda data, qty, logNpmin : (data["log_"+qty]>data["logMpart"]+logNpmin)

zSelection = lambda data, zmin, zmax : (data["redshift"]>zmin)&(data["redshift"]<zmax)

nSelection = lambda data, NminCount, cos : (data['dN_counts_'+cos]>NminCount)

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
zSel = zSelection( data, zmin, zmax )
# mass selection
mSel = mSelection(data, qty, logNpmin)
mSel2_inter = (data["log_mvir"]<13.2) & (data["redshift"]>0.)
mSel2 = (mSel2_inter==False)
# minimum number counts selection
nSelCen = nSelection(data, NminCount, cos )
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


bias_all = n.array([
f_BH(sigma, 0.3241+0.0004, 0.897+0.006, 0.624+0.025, 1.589+0.03), 
f_BH(sigma, 0.3241+0.0004, 0.897+0.006, 0.624+0.025, 1.589-0.03), 
f_BH(sigma, 0.3241+0.0004, 0.897-0.006, 0.624+0.025, 1.589-0.03), 
f_BH(sigma, 0.3241+0.0004, 0.897-0.006, 0.624+0.025, 1.589+0.03), 
f_BH(sigma, 0.3241-0.0004, 0.897-0.006, 0.624+0.025, 1.589-0.03), 
f_BH(sigma, 0.3241-0.0004, 0.897-0.006, 0.624+0.025, 1.589+0.03), 
f_BH(sigma, 0.3241-0.0004, 0.897+0.006, 0.624+0.025, 1.589-0.03),
f_BH(sigma, 0.3241-0.0004, 0.897+0.006, 0.624+0.025, 1.589+0.03),
f_BH(sigma, 0.3241+0.0004, 0.897+0.006, 0.624-0.025, 1.589+0.03), 
f_BH(sigma, 0.3241+0.0004, 0.897+0.006, 0.624-0.025, 1.589-0.03), 
f_BH(sigma, 0.3241+0.0004, 0.897-0.006, 0.624-0.025, 1.589-0.03), 
f_BH(sigma, 0.3241+0.0004, 0.897-0.006, 0.624-0.025, 1.589+0.03), 
f_BH(sigma, 0.3241-0.0004, 0.897-0.006, 0.624-0.025, 1.589-0.03), 
f_BH(sigma, 0.3241-0.0004, 0.897-0.006, 0.624-0.025, 1.589+0.03), 
f_BH(sigma, 0.3241-0.0004, 0.897+0.006, 0.624-0.025, 1.589-0.03),
f_BH(sigma, 0.3241-0.0004, 0.897+0.006, 0.624-0.025, 1.589+0.03)
])
bi_max = n.log10(n.max(bias_all,axis = 0))
bi_min = n.log10(n.min(bias_all,axis = 0))

bias_all = n.array([
f_ST(sigma, 0.04235+0.00037, 1.70219+0.01035, 0.83118+0.03672), 
f_ST(sigma, 0.04235+0.00037, 1.70219+0.01035, 0.83118-0.03672), 
f_ST(sigma, 0.04235+0.00037, 1.70219-0.01035, 0.83118-0.03672), 
f_ST(sigma, 0.04235+0.00037, 1.70219-0.01035, 0.83118+0.03672), 
f_ST(sigma, 0.04235-0.00037, 1.70219-0.01035, 0.83118-0.03672), 
f_ST(sigma, 0.04235-0.00037, 1.70219-0.01035, 0.83118+0.03672), 
f_ST(sigma, 0.04235-0.00037, 1.70219+0.01035, 0.83118-0.03672),
f_ST(sigma, 0.04235-0.00037, 1.70219+0.01035, 0.83118+0.03672)
])
bi2_max = n.log10(n.max(bias_all,axis = 0))
bi2_min =n.log10( n.min(bias_all,axis = 0))



p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
# log_mvir, logsigM1, logNu, log_MF, log_MF_c, redshift  = log_mvir[ok], logsig[ok], lognu[ok], log_MF[ok], log_MF_c[ok], data['redshift'][ok]

ok = (zSel) & (mSel) & (mSel2) & (nSelCen)&(data["boxName"]!='DS_8Gpc')
print "distinct"
print "N distinct z=0 ",  n.sum(data['dN_counts_cen'][ok])
print "M range ",  n.min(data['log_mvir_min'][ok]),  n.max(data['log_mvir_max'][ok])
x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_cen"][ok]**2. + data["dN_counts_cen"][ok]**(-1.))**(0.5)
p.errorbar(x_data, y_data, yerr = y_data_err, rasterized=True, fmt='none', label='distinct halos z=0', lw=1)
p.fill_between(-n.log10(sigma), y1=bi_min, y2=bi_max, color='k',alpha=0.3)

cos = 'sat'
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )
nSelSat = nSelection(data, NminCount, cos )
ok = (zSel) & (mSel)& (mSel2)  & (nSelSat)&(data["boxName"]!='DS_8Gpc')
print "sat"
print "N distinct z=0 ",  n.sum(data['dN_counts_sat'][ok])
print "M range ",  n.min(data['log_mvir_min'][ok]),  n.max(data['log_mvir_max'][ok])
x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_"+cos][ok]**2. + data["dN_counts_"+cos][ok]**(-1.))**(0.5)
p.errorbar(x_data, y_data, yerr = y_data_err, rasterized=True, fmt='none', label='satellite subhalos z=0', lw=1)
sigs = n.arange(-0.5,.6, 0.01)
p.fill_between(-n.log10(sigma), y1=bi2_min, y2=bi2_max, color='k',alpha=0.3, label='model')

p.xlabel(r'$\log_{10}(\sigma^{-1})$')
p.ylabel(r'$\log_{10}\left[ \frac{M}{\rho_m} \frac{dn}{d\ln M} \left|\frac{d\ln M }{d\ln \sigma}\right|\right] $') 
 # log$_{10}[ n(>M)]')
gl = p.legend(loc=0, fontsize=12)
gl.set_frame_on(False)
p.ylim((-3., -0.4))
p.xlim((-0.5, 0.5))
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],"mvir-AL-z0-differential-function-data-xSigma.png"))
p.clf()

#=======================
#=======================
#=======================
#=======================
# RESIDUAL PLOT
#=======================
#=======================
#=======================
#=======================


#=======================
#=======================
cos = 'cen'
#=======================
#=======================

x_model = -n.log10(sigma)
y_model = f_BH(sigma, 0.3241, 0.897, 0.624, 1.589)
model_interpol = interp1d(x_model, y_model)
y_model_min = 10**bi_min
y_model_max = 10**bi_max

# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )

# redshift selection
zSel = zSelection( data, zmin, zmax )
# mass selection
mSel = mSelection(data, qty, logNpmin)
# minimum number counts selection
nSelCen = nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (mSel2) & (nSelCen)
# selection per box :
ok = (zSel) & (mSel) & (mSel2) & (nSelCen)&(data["boxName"]!='DS_8Gpc')

x_data = logsig
y_data = log_MF
y_data_err = (data["std90_pc_cen"]**2. + data["dN_counts_cen"]**(-1.))**(0.5)

MD_sel_fun=lambda name : (ok)&(data["boxName"]==name)
MDnames= n.array(['M04', 'M10', 'M25','M40','M25n','M40n'])
MDconv= n.array(['MD_0.4Gpc', 'MD_1Gpc', 'MD_2.5Gpc','MD_4Gpc','MD_2.5GpcNW','MD_4GpcNW'])
MDsels=n.array([MD_sel_fun(name) for name in MDconv])

f_diff_fun = lambda MDs:  10**y_data[MDs] / model_interpol(x_data[MDs])
f_diffs = n.array([f_diff_fun(MD) for MD in MDsels])

print "================================"

# now the plots
p.figure(0,(6,6))
p.axes([0.15,0.12,0.75,0.75])
for index, fd in enumerate(f_diffs):
	inTol = (abs(10**fd-1)<tolerance)
	print index
	if len(fd)>0:
		p.errorbar(x_data[MDsels[index]], fd, yerr = y_data_err[MDsels[index]] , rasterized=True, fmt='none', label=MDnames[index])
		print len(inTol.nonzero()[0]), len(fd), 100.*len(inTol.nonzero()[0])/ len(fd)

p.fill_between(x_model, y1=y_model_min/y_model, y2=y_model_max/y_model, color='k',alpha=0.2, label='model')
#p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
#p.axhline(0.99,c='k',ls='--')
p.xlabel(r'$\log_{10}(\sigma^{-1})$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
p.xlim((-0.5, 0.5))
p.ylim((0.9,1.1))
#p.yscale('log')
p.grid()
p.title('distinct')
p.savefig(join(dir,"fit-"+cos+"-differential-function-residual-log.png"))
p.clf()

print "DISTINCT mean and std of the residuals", n.mean(n.hstack((f_diffs)) - 1.), n.std(n.hstack((f_diffs)) - 1.)


#=======================
#=======================
cos = 'sat'
#=======================
#=======================

x_model = -n.log10(sigma)
y_model = f_ST(sigma, 0.04235, 1.70219, 0.83118) 
model_interpol = interp1d(x_model, y_model)
y_model_min = 10**bi2_min
y_model_max = 10**bi2_max

# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )

# redshift selection
zSel = zSelection( data, zmin, zmax )
# mass selection
mSel = mSelection(data, qty, logNpmin)
# minimum number counts selection
nSelCen = nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (mSel2) & (nSelCen)
# selection per box :
ok = (zSel) & (mSel) & (mSel2) & (nSelCen)&(data["boxName"]!='DS_8Gpc')

x_data = logsig
y_data = log_MF
y_data_err = (data["std90_pc_cen"]**2. + data["dN_counts_cen"]**(-1.))**(0.5)

MD_sel_fun=lambda name : (ok)&(data["boxName"]==name)
MDnames= n.array(['M04', 'M10', 'M25','M40','M25n','M40n'])
MDconv= n.array(['MD_0.4Gpc', 'MD_1Gpc', 'MD_2.5Gpc','MD_4Gpc','MD_2.5GpcNW','MD_4GpcNW'])
MDsels=n.array([MD_sel_fun(name) for name in MDconv])

f_diff_fun = lambda MDs:  10**y_data[MDs] / model_interpol(x_data[MDs])
f_diffs = n.array([f_diff_fun(MD) for MD in MDsels])

print "================================"

# now the plots
p.figure(0,(6,6))
p.axes([0.15,0.12,0.75,0.75])
for index, fd in enumerate(f_diffs):
	inTol = (abs(10**fd-1)<tolerance)
	print index
	if len(fd)>0:
		p.errorbar(x_data[MDsels[index]], fd, yerr = y_data_err[MDsels[index]] , rasterized=True, fmt='none', label=MDnames[index])
		print len(inTol.nonzero()[0]), len(fd), 100.*len(inTol.nonzero()[0])/ len(fd)

p.fill_between(x_model, y1=y_model_min/y_model, y2=y_model_max/y_model, color='k',alpha=0.2, label='model')
#p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
#p.axhline(0.99,c='k',ls='--')
p.xlabel(r'$\log_{10}(\sigma^{-1})$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
p.xlim((-0.5, 0.5))
p.ylim((0.9,1.1))
#p.yscale('log')
p.grid()
p.title('satelitte')
p.savefig(join(dir,"fit-"+cos+"-differential-function-residual-log.png"))
p.clf()

print "SUB mean and std of the residuals", n.mean(n.hstack((f_diffs)) - 1.), n.std(n.hstack((f_diffs)) - 1.)

