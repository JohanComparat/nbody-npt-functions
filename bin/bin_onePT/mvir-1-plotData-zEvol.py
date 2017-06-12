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

from scipy.interpolate import interp1d
#Quantity studied
qty = "mvir"
# working directory
dir = join(os.environ['MVIR_DIR'])
# loads summary file
data = fits.open( join(dir, qty+"_summary.fits"))[1].data

NminCount = 1000
logNpmin = 3

zmin = -0.01
zmax = 2.3

tolerance = 0.03

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

bias_all = n.array([
lib.f_BH(lib.hmf.sigma, 0.28074+0.00151, 0.90343+0.00724,  0.64031+0.02639, 1.69561+0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074+0.00151, 0.90343+0.00724, 0.64031+0.02639, 1.69561-0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074+0.00151, 0.90343-0.00724, 0.64031+0.02639, 1.69561-0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074+0.00151, 0.90343-0.00724, 0.64031+0.02639, 1.69561+0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074-0.00151, 0.90343-0.00724, 0.64031+0.02639, 1.69561-0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074-0.00151, 0.90343-0.00724, 0.64031+0.02639, 1.69561+0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074-0.00151, 0.90343+0.00724, 0.64031+0.02639, 1.69561-0.03826),
lib.f_BH(lib.hmf.sigma, 0.28074-0.00151, 0.90343+0.00724, 0.64031+0.02639, 1.69561+0.03826),
lib.f_BH(lib.hmf.sigma, 0.28074+0.00151, 0.90343+0.00724,  0.64031-0.02639, 1.69561+0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074+0.00151, 0.90343+0.00724, 0.64031-0.02639, 1.69561-0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074+0.00151, 0.90343-0.00724, 0.64031-0.02639, 1.69561-0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074+0.00151, 0.90343-0.00724, 0.64031-0.02639, 1.69561+0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074-0.00151, 0.90343-0.00724, 0.64031-0.02639, 1.69561-0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074-0.00151, 0.90343-0.00724, 0.64031-0.02639, 1.69561+0.03826), 
lib.f_BH(lib.hmf.sigma, 0.28074-0.00151, 0.90343+0.00724, 0.64031-0.02639, 1.69561-0.03826),
lib.f_BH(lib.hmf.sigma, 0.28074-0.00151, 0.90343+0.00724, 0.64031-0.02639, 1.69561+0.03826)
])
bi_max = n.log10(n.max(bias_all,axis = 0))
bi_min = n.log10(n.min(bias_all,axis = 0))

bias_all = n.array([
lib.f_ST(lib.hmf.sigma, 0.04235+0.00037, 1.70219+0.01035, 0.83118+0.03672), 
lib.f_ST(lib.hmf.sigma, 0.04235+0.00037, 1.70219+0.01035, 0.83118-0.03672), 
lib.f_ST(lib.hmf.sigma, 0.04235+0.00037, 1.70219-0.01035, 0.83118-0.03672), 
lib.f_ST(lib.hmf.sigma, 0.04235+0.00037, 1.70219-0.01035, 0.83118+0.03672), 
lib.f_ST(lib.hmf.sigma, 0.04235-0.00037, 1.70219-0.01035, 0.83118-0.03672), 
lib.f_ST(lib.hmf.sigma, 0.04235-0.00037, 1.70219-0.01035, 0.83118+0.03672), 
lib.f_ST(lib.hmf.sigma, 0.04235-0.00037, 1.70219+0.01035, 0.83118-0.03672),
lib.f_ST(lib.hmf.sigma, 0.04235-0.00037, 1.70219+0.01035, 0.83118+0.03672)
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
#p.plot(lib.X, n.log10(lib.ftC16), 'k--', label='fit', lw=2)
p.fill_between(-n.log10(lib.hmf.sigma), y1=bi_min, y2=bi_max, color='k',alpha=0.3)

cos = 'sat'
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )
nSelSat = lib.nSelection(data, NminCount, cos )
ok = (zSel) & (mSel)& (mSel2)  & (nSelSat)&(data["boxName"]!='DS_8Gpc')
print "sat"
print "N distinct z<2.3 ",  n.sum(data['dN_counts_sat'][ok])
print "M range ",  n.min(data['log_mvir_min'][ok]),  n.max(data['log_mvir_max'][ok])
x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_"+cos][ok]**2. + data["dN_counts_"+cos][ok]**(-1.))**(0.5)
p.errorbar(x_data, y_data, yerr = y_data_err, rasterized=True, fmt='none', label='satellite subhalos z<2.3', lw=1)
sigs = n.arange(-0.5,.6, 0.01)
#p.plot(lib.X, n.log10(lib.ftC16st_sat), 'k--', lw=2)
p.fill_between(-n.log10(lib.hmf.sigma), y1=bi2_min, y2=bi2_max, color='k',alpha=0.3, label='model')

p.xlabel(r'$\log_{10}(\sigma^{-1})$')
p.ylabel(r'$\log_{10}\left[ \frac{M}{\rho_m} \frac{dn}{d\ln M} \left|\frac{d\ln M }{d\ln \sigma}\right|\right] $') 
 # log$_{10}[ n(>M)]')
gl = p.legend(loc=0, fontsize=12)
gl.set_frame_on(False)
p.ylim((-3., -0.4))
p.xlim((-0.5, 0.5))
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],"mvir-AL-z_lt_23-differential-function-data-xSigma.png"))
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

x_model = -n.log10(lib.hmf.sigma)
y_model = lib.f_BH(lib.hmf.sigma, 0.28074, 0.90343, 0.64031, 1.695616)
model_interpol = interp1d(x_model, y_model)
y_model_min = 10**bi_min
y_model_max = 10**bi_max

# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )

# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, logNpmin)
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
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
p.axes([0.17,0.17,0.75,0.75])
for index, fd in enumerate(f_diffs):
	inTol = (abs(10**fd-1)<tolerance)
	print index
	if len(fd)>0:
		p.scatter(x_data[MDsels[index]], fd, c= data["redshift"][MDsels[index]], s=20)
		#p.errorbar(x_data[MDsels[index]], fd, yerr = y_data_err[MDsels[index]] , rasterized=True, fmt='none', label=MDnames[index])
		#print len(inTol.nonzero()[0]), len(fd), 100.*len(inTol.nonzero()[0])/ len(fd)

p.fill_between(x_model, y1=y_model_min/y_model, y2=y_model_max/y_model, color='k',alpha=0.2, label='model')
#p.axhline(1.01,c='k',ls='--',label=r'syst $\pm1\%$')
#p.axhline(0.99,c='k',ls='--')
p.xlabel(r'$\log_{10}(\sigma^{-1})$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
p.colorbar()
gl.set_frame_on(False)
p.xlim((-0.5, 0.5))
p.ylim((0.1,1.9))
#p.yscale('log')
p.grid()
p.savefig(join(dir,"fit-"+cos+"-differential-function-residual-log-zlt23.png"))
p.clf()

print "DISTINCT mean and std of the residuals", n.mean(n.hstack((f_diffs)) - 1.), n.std(n.hstack((f_diffs)) - 1.)

sys.exit()
#=======================
#=======================
cos = 'sat'
#=======================
#=======================

x_model = -n.log10(lib.hmf.sigma)
y_model = lib.f_ST(lib.hmf.sigma, 0.04235, 1.70219, 0.83118) 
model_interpol = interp1d(x_model, y_model)
y_model_min = 10**bi2_min
y_model_max = 10**bi2_max

# y coordinates
ff = mvir *  data["dNdlnM_"+cos] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
ff_c = mvir *  data["dNdlnM_"+cos+"_c"] / data["rhom"]  / abs(data["dlnsigmaMdlnM"]) 
log_MF = n.log10( ff )
log_MF_c = n.log10(  ff_c )

# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, logNpmin)
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
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
p.axes([0.17,0.17,0.75,0.75])
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
p.savefig(join(dir,"fit-"+cos+"-differential-function-residual-log.png"))
p.clf()

print "SUB mean and std of the residuals", n.mean(n.hstack((f_diffs)) - 1.), n.std(n.hstack((f_diffs)) - 1.)


sys.exit()


x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_"+cos][ok]**2. + data["dN_counts_"+cos][ok]**(-1.))**(0.5)
p.errorbar(x_data, y_data, yerr = y_data_err, rasterized=True, fmt='none', label='satellite subhalos z=0', lw=1)
sigs = n.arange(-0.5,.6, 0.01)
#p.plot(lib.X, n.log10(lib.ftC16st_sat), 'k--', lw=2)
p.fill_between(-n.log10(lib.hmf.sigma), y1=bi2_min, y2=bi2_max, color='k',alpha=0.3, label='model')

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

