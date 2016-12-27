import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
import lib_functions_1pt as lib
import os
import sys

# data modules
import glob
import sys
import astropy.io.fits as fits
import os
from os.path import join
import cPickle

# numerical modules
import numpy as n
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.interpolate import interp2d
from scipy.stats import norm
from scipy.interpolate import griddata

# plotting modules
import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

# mass function theory
from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

#lib.covariance_factor
#lib.f_BH(sigma, 0.333, 0.788, 0.807, 1.795)
A0=0.333
a0=0.786
p0=0.807
q0=1.795
bias = lambda sigma, a0, p0, q0 : lib.b_BH(sigma, a0, p0, q0)
fsigma = lambda sigma : lib.f_BH(sigma, A0, a0, p0, q0)

# diagonal error
dn_L04 = lambda sigma, a, p, q :  (((bias(sigma, a, p, q) * lib.hmf.growth_factor)**2. * (lib.covariance_factor[0]))**2. + lib.shot_noise(sigma, 400.**3.)  )**0.5
dn_L10 = lambda sigma, a, p, q :  (((bias(sigma, a, p, q) * lib.hmf.growth_factor)**2. * (lib.covariance_factor[1]) )**2. + lib.shot_noise(sigma, 1000.**3.) )**0.5  
dn_L25 = lambda sigma, a, p, q : (((bias(sigma, a, p, q) * lib.hmf.growth_factor)**2. * (lib.covariance_factor[2]) )**2. + lib.shot_noise(sigma, 2500.**3.) )**0.5  
dn_L40 = lambda sigma, a, p, q :  (((bias(sigma, a, p, q) * lib.hmf.growth_factor)**2. * (lib.covariance_factor[3]) )**2. + lib.shot_noise(sigma, 4000.**3.) )**0.5  

# off diagonal error 
dn_cov_L04 = lambda s1, s2, a, p, q : (dn_L04(s1, a, p, q)*dn_L04(s2, a, p, q))**0.5 
dn_cov_L10 = lambda s1, s2, a, p, q : (dn_L10(s1, a, p, q)*dn_L10(s2, a, p, q) )**0.5
dn_cov_L25 = lambda s1, s2, a, p, q : (dn_L25(s1, a, p, q)*dn_L25(s2, a, p, q) )**0.5
dn_cov_L40 = lambda s1, s2, a, p, q : (dn_L40(s1, a, p, q)*dn_L40(s2, a, p, q) )**0.5

# opens the data 

#Quantity studied
qty = "mvir"
# working directory
dir = join(os.environ['MVIR_DIR'])
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
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')

x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_cen"][ok]**2. + data["dN_counts_cen"][ok]**(-1.))**(0.5)

ps = n.array([0.333, 0.794, 0.247])
log_fsigma = lambda logsigma, A, a, p : n.log10(lib.f_BH(10**-logsigma, A, a, p, 1.))
print "ST01 fit ----------------------"
pOpt, pCov=curve_fit(log_fsigma, x_data, y_data, ps, sigma = y_data_err, maxfev=50000000)#, bounds=boundaries)
chi2 = n.sum(((log_fsigma(x_data, pOpt[0], pOpt[1], pOpt[2])-y_data)/y_data_err)**2. ) 
ndof = (len(x_data) - len(ps)) 
print "best params=", pOpt
print "err=", pCov.diagonal()**0.5
print "chi2 ", chi2, ndof, chi2/ndof
print "---------------------------------------------------"
pOpt_ST01 = pOpt
pErr_ST01 = pCov.diagonal()**0.5

n.savetxt(join(os.environ['MVIR_DIR'],"mvirFunction_parameters_ST01_MFonly_fit.txt"), n.transpose([pOpt_ST01, pErr_ST01]), header="A a p")

print "BATT 2011"
ps = n.array([A0, a0, p0, q0])
log_fsigma = lambda logsigma, A, a, p, q : n.log10(lib.f_BH(10**-logsigma, A, a, p, q))

pOpt, pCov=curve_fit(log_fsigma, x_data, y_data, ps, sigma = y_data_err, maxfev=50000000)#, bounds=boundaries)
chi2 = n.sum(((log_fsigma(x_data, pOpt[0], pOpt[1], pOpt[2], pOpt[3])-y_data)/y_data_err)**2. ) 
ndof = (len(x_data) - len(ps)) 
print "best params=", pOpt
print "err=", pCov.diagonal()**0.5
print "chi2 ", chi2, ndof, chi2/ndof
print "---------------------------------------------------"
A1, a1, p1, q1 = pOpt
A1_err, a1_err, p1_err, q1_err = pCov.diagonal()**0.5

n.savetxt(join(os.environ['MVIR_DIR'],"mvirFunction_parameters_BA11_MFonly_fit.txt"), n.transpose([pOpt, pCov.diagonal()**0.5]), header="A a p q")

def plotSel(MDsel, label):
	x_data = logsig[ok & MDsel]
	y_data = log_MF[ok & MDsel]
	y_data_err = (data["std90_pc_cen"][ok & MDsel]**2. + data["dN_counts_cen"][ok & MDsel]**(-1.))**(0.5)
	y_model = log_fsigma(x_data, pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	if len(x_data)>0:
		p.errorbar(x_data, 10**(y_data-y_model), yerr = y_data_err , rasterized=True, fmt='none', label=label)


p.figure(0,(6,3))
p.axes([0.17,0.17,0.75,0.6])
plotSel(MD04, "MD04")
plotSel(MD10, "MD10")
plotSel(MD25, "MD25")
plotSel(MD40, "MD40")
plotSel(MD25NW, "MD25NW")
plotSel(MD40NW, "MD40NW")
p.axhline(1.025,c='k',ls='--',label=r'$\pm2.5\%$')
p.axhline(0.975,c='k',ls='--')
p.xlabel(r'$log_{10}(\sigma^{-1})$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
#p.xlim((-0.7,0.6))
p.ylim((0.9,1.15))
#p.yscale('log')
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],"fit-BA11-"+cos+"-differential-function-residual.png"))
p.clf()

ok2 = (zSel) & (nSelCen)
def plotAll(MDsel, label):
	x_data = logsig[ok2 & MDsel]
	x_plot = data["log_mvir"][ok2 & MDsel]
	y_data = log_MF[ok2 & MDsel]
	y_data_err = (data["std90_pc_cen"][ok2 & MDsel]**2. + data["dN_counts_cen"][ok2 & MDsel]**(-1.))**(0.5)
	y_model = log_fsigma(x_data, pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	if len(x_data)>0:
		p.errorbar(x_plot, 10**(y_data-y_model), yerr = y_data_err , rasterized=True, label=label)

		
p.figure(0,(6,3))
p.axes([0.17,0.17,0.75,0.6])
plotAll(MD04, "MD04")
plotAll(MD10, "MD10")
plotAll(MD25, "MD25")
plotAll(MD25NW, "MD25NW")
plotAll(MD40, "MD40")
plotAll(MD40NW, "MD40NW")

p.axhline(0.8,c='k',ls='--')
p.axhline(0.9,c='k',ls='--')
p.axhline(0.95,c='k',ls='--')
p.axhline(0.97,c='k',ls='--')
p.xlabel(r'$log_{10}(M_h)$')
p.ylabel(r'data/model') 
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
p.xlim((9.5, 16))
p.ylim((0.75,1.05))
#p.yscale('log')
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],"fit-BA11-"+cos+"-differential-function-residual-incompleteness.png"))
p.clf()

fileC = n.array(glob.glob( join(os.environ['MD_DIR'],"MD_*Gpc*", "properties", qty,"out_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MD_DIR'],"MD_*Gpc*","properties", qty,"out_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MD_DIR'],"MD_*Gpc*","properties", qty,"out_*_Satellite_JKresampling.pkl")))

# redshift 0 data
iis = [8, 13, 29, 31, 60, -10]#[-1, -2, -4, -9, -22, 3]

def get_COV(fileCov, binFile):
	boxZN = float(os.path.basename(fileCov).split('_')[1])
	print boxZN
	hf, boxLength, boxName, boxRedshift, logmp, boxLengthComoving, massCorrection = lib.get_basic_info(fileCov, boxZN, delta_wrt='mean')
	bins = n.log10( 10**n.loadtxt(binFile) * massCorrection )
	logmass = ( bins[1:]  + bins[:-1] )/2.
	mass = 10**logmass
	dX = ( 10**bins[1:]  - 10**bins[:-1] )
	dlnbin = (bins[1:]  - bins[:-1])*n.log(10)
	m2sigma = interp1d(hf.M, hf.sigma )
	sigma_i = m2sigma( mass )
	data_i=cPickle.load(open(fileCov,'r'))
	counts_i = n.sum(data_i, axis=0)
	ok = (counts_i>10)
	data = data_i.T[ok].T
	counts = n.sum(data, axis=0)
	sigma = sigma_i[ok]
	count_matrix = n.outer(counts, counts)/1000.
	cv = (n.cov(data.T, ddof=0)/count_matrix)**0.5
	ctotal = cv + count_matrix**(-0.5)
	xcv, ycv = n.meshgrid(sigma, sigma)
	s1_i = n.hstack((n.log10(xcv)))
	s2_i = n.hstack((n.log10(ycv)))
	val_i = n.hstack((n.log10(ctotal)))
	ok = (val_i != n.inf) & (n.isnan(val_i) == False )
	return s1_i[ok], s2_i[ok], val_i[ok]


model = {"MD_0.4Gpc": dn_cov_L04, "MD_1Gpc": dn_cov_L10, "MD_2.5Gpc": dn_cov_L25, "MD_4Gpc": dn_cov_L40,"MD_2.5GpcNW": dn_cov_L25, "MD_4GpcNW": dn_cov_L40 }

qs = n.arange(0,4,0.01)

ps0 = [ q1]
error =  0.02 # dex 
log_cov_04_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L04( 10**logs1, 10**logs2, a1, ps[0], q1 ) )
log_cov_10_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L10( 10**logs1, 10**logs2, a1, ps[0], q1 ) )
log_cov_25_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L25( 10**logs1, 10**logs2, a1, ps[0], q1 ) )
log_cov_40_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L40( 10**logs1, 10**logs2, a1, ps[0], q1 ) )

# MD 04
index=0
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_04_ps(log_s1_data, log_s2_data, ps) - cov_data)**2.  /( error ))/(len(cov_data) - len(ps0))
chi2fun_pp = lambda ps :  (log_cov_04_ps(log_s1_data, log_s2_data, ps) - cov_data)**2.  /( error )/(len(cov_data) - len(ps0))

res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_04 = res.x
pCov_04 = res.direc
print "04, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5

print "====================="
chi2s_qs_04 = n.array([chi2fun([qq]) for qq in qs ])



# MD 10
index=1
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_10_ps(log_s1_data, log_s2_data, ps) - cov_data)**2. /( error ))/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_10 = res.x
pCov_10 = res.direc
print "10, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="
chi2s_qs_10 = n.array([chi2fun([qq]) for qq in qs ])

# MD 25
index=2
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_25_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_25 = res.x
pCov_25 = res.direc
print "25, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="
chi2s_qs_25 = n.array([chi2fun([qq]) for qq in qs ])

# MD 25 NW
index=3
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_25_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_25NW = res.x
pCov_25NW = res.direc
print "25NW, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="
chi2s_qs_25NW = n.array([chi2fun([qq]) for qq in qs ])

# MD 40
index=4
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_40_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_40 = res.x
pCov_40 = res.direc
print "40, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="
chi2s_qs_40 = n.array([chi2fun([qq]) for qq in qs ])

# MD 40
index=5
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_40_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_40NW = res.x
pCov_40NW = res.direc
print "40NW, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="

chi2s_qs_40NW = n.array([chi2fun([qq]) for qq in qs ])

pOpts = n.array([pOpt_04, pOpt_10, pOpt_25, pOpt_25NW, pOpt_40, pOpt_40NW ])
pErrs = n.array([pCov_04.diagonal()**0.5, pCov_10.diagonal()**0.5, pCov_25.diagonal()**0.5, pCov_25NW.diagonal()**0.5, pCov_40.diagonal()**0.5, pCov_40NW.diagonal()**0.5 ])
print "best params=", pOpts.T
print "err=", pErrs.T

p.figure(8)
p.plot(qs, chi2s_qs_04, label='04')
p.plot(qs, chi2s_qs_10, label='10')
p.plot(qs, chi2s_qs_25, label='25')
p.plot(qs, chi2s_qs_25NW, label='25')
p.plot(qs, chi2s_qs_40, label='40')
p.plot(qs, chi2s_qs_40NW, label='40NW')
p.axvline(p1+p1_err, c='k', ls='dashed', label='MF fit')
p.axvline(p1-p1_err, c='k', ls='dashed')
p.legend(loc=0, frameon=False)
p.xlim((0,3))
p.ylim((0,5))
p.xlabel('p')
p.ylabel(r'$\chi^2/ndof$')
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],"covariance","p_parameter_constrain.png"))
p.clf()


qs = n.arange(0,4,0.01)

ps0 = [ q1]
error =  0.02 # dex 
log_cov_04_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L04( 10**logs1, 10**logs2, a1, p1, ps[0] ) )
log_cov_10_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L10( 10**logs1, 10**logs2, a1, p1, ps[0] ) )
log_cov_25_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L25( 10**logs1, 10**logs2, a1, p1, ps[0] ) )
log_cov_40_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L40( 10**logs1, 10**logs2, a1, p1, ps[0] ) )

# MD 04
index=0
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_04_ps(log_s1_data, log_s2_data, ps) - cov_data)**2.  /( error ))/(len(cov_data) - len(ps0))
chi2fun_pp = lambda ps :  (log_cov_04_ps(log_s1_data, log_s2_data, ps) - cov_data)**2.  /( error )/(len(cov_data) - len(ps0))

res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_04 = res.x
pCov_04 = res.direc
print "04, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="
chi2s_qs_04 = n.array([chi2fun([qq]) for qq in qs ])

# MD 10
index=1
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_10_ps(log_s1_data, log_s2_data, ps) - cov_data)**2. /( error ))/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_10 = res.x
pCov_10 = res.direc
print "10, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="
chi2s_qs_10 = n.array([chi2fun([qq]) for qq in qs ])

# MD 25
index=2
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_25_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_25 = res.x
pCov_25 = res.direc
print "25, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="
chi2s_qs_25 = n.array([chi2fun([qq]) for qq in qs ])

# MD 25 NW
index=3
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_25_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_25NW = res.x
pCov_25NW = res.direc
print "25NW, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="
chi2s_qs_25NW = n.array([chi2fun([qq]) for qq in qs ])

# MD 40
index=4
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_40_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_40 = res.x
pCov_40 = res.direc
print "40, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="
chi2s_qs_40 = n.array([chi2fun([qq]) for qq in qs ])

# MD 40
index=5
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_40_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_40NW = res.x
pCov_40NW = res.direc
print "40NW, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "====================="

chi2s_qs_40NW = n.array([chi2fun([qq]) for qq in qs ])

pOpts = n.array([pOpt_04, pOpt_10, pOpt_25, pOpt_25NW, pOpt_40, pOpt_40NW ])
pErrs = n.array([pCov_04.diagonal()**0.5, pCov_10.diagonal()**0.5, pCov_25.diagonal()**0.5, pCov_25NW.diagonal()**0.5, pCov_40.diagonal()**0.5, pCov_40NW.diagonal()**0.5 ])
print "best params=", pOpts.T
print "err=", pErrs.T

p.figure(8)
p.plot(qs, chi2s_qs_04, label='04')
p.plot(qs, chi2s_qs_10, label='10')
p.plot(qs, chi2s_qs_25, label='25')
p.plot(qs, chi2s_qs_25NW, label='25')
p.plot(qs, chi2s_qs_40, label='40')
p.plot(qs, chi2s_qs_40NW, label='40NW')
p.axvline(q1+q1_err, c='k', ls='dashed', label='MF fit')
p.axvline(q1-q1_err, c='k', ls='dashed')
p.legend(loc=0, frameon=False)
p.xlim((0,3))
p.ylim((0,5))
p.xlabel('q')
p.ylabel(r'$\chi^2/ndof$')
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],"covariance","q_parameter_constrain.png"))
p.clf()


print "------------------------------------------"
print "------------------------------------------"
print "------------------------------------------"
print "------2d 2d 2d 2d 2d 2d 2d2d 2d  -----"
print "------------------------------------------"
print "------------------------------------------"
print "------------------------------------------"


ps_i, qs_i = n.meshgrid(n.arange(0,4,0.01),n.arange(0,4,0.01))
ps = n.hstack((ps_i))
qs = n.hstack((qs_i))

ps0 = [ p1, q1]
error =  0.02 # dex 
log_cov_04_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L04( 10**logs1, 10**logs2, a1, ps[0], ps[1] ) )
log_cov_10_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L10( 10**logs1, 10**logs2, a1, ps[0], ps[1] ) )
log_cov_25_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L25( 10**logs1, 10**logs2, a1, ps[0], ps[1] ) )
log_cov_40_ps = lambda logs1, logs2, ps : n.log10( dn_cov_L40( 10**logs1, 10**logs2, a1, ps[0], ps[1] ) )

# MD 04
index=0
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_04_ps(log_s1_data, log_s2_data, ps) - cov_data)**2.  /( error ))/(len(cov_data) - len(ps0))
chi2fun_pp = lambda ps :  (log_cov_04_ps(log_s1_data, log_s2_data, ps) - cov_data)**2.  /( error )/(len(cov_data) - len(ps0))

res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_04 = res.x
pCov_04 = res.direc
print "04, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "chi2 ", chi2fun(res.x)*(len(cov_data) - len(ps0)), len(cov_data) - len(ps0), chi2fun(res.x)

print "====================="
#chi2s_qs_04 = n.array([n.array([chi2fun([pp, qq]) for qq in qs ]) for pp in ps ])

# MD 10
index=1
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_10_ps(log_s1_data, log_s2_data, ps) - cov_data)**2. /( error ))/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_10 = res.x
pCov_10 = res.direc
print "10, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "chi2 ", chi2fun(res.x)*(len(cov_data) - len(ps0)), len(cov_data) - len(ps0), chi2fun(res.x)

print "====================="
#chi2s_qs_10 = n.array([n.array([chi2fun([pp, qq]) for qq in qs ]) for pp in ps ])
# MD 25
index=2
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_25_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_25 = res.x
pCov_25 = res.direc
print "25, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "chi2 ", chi2fun(res.x)*(len(cov_data) - len(ps0)), len(cov_data) - len(ps0), chi2fun(res.x)

print "====================="
#chi2s_qs_25 = n.array([n.array([chi2fun([pp, qq]) for qq in qs ]) for pp in ps ])
# MD 25 NW
index=3
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_25_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_25NW = res.x
pCov_25NW = res.direc
print "25NW, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "chi2 ", chi2fun(res.x)*(len(cov_data) - len(ps0)), len(cov_data) - len(ps0), chi2fun(res.x)

print "====================="
#chi2s_qs_25NW = n.array([n.array([chi2fun([pp, qq]) for qq in qs ]) for pp in ps ])
# MD 40
index=4
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_40_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_40 = res.x
pCov_40 = res.direc
print "40, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "chi2 ", chi2fun(res.x)*(len(cov_data) - len(ps0)), len(cov_data) - len(ps0), chi2fun(res.x)

print "====================="
#chi2s_qs_40 = n.array([n.array([chi2fun([pp, qq]) for qq in qs ]) for pp in ps ])
# MD 40
index=5
log_s1_data, log_s2_data, cov_data = get_COV(fileC[iis[index]], fileB[iis[index]])
chi2fun = lambda ps : n.sum( (log_cov_40_ps(log_s1_data, log_s2_data, ps) - cov_data)**2./( error ) )/(len(cov_data) - len(ps0))
res = minimize(chi2fun, ps0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
pOpt_40NW = res.x
pCov_40NW = res.direc
print "40NW, init chi2=", chi2fun(ps0)
print "best params=", res.x
print "err=", res.direc.diagonal()**0.5
print "chi2 ", chi2fun(res.x)*(len(cov_data) - len(ps0)), len(cov_data) - len(ps0), chi2fun(res.x)

print "====================="

#chi2s_qs_40NW = n.array([n.array([chi2fun([pp, qq]) for qq in qs ]) for pp in ps ])

pOpts = n.array([pOpt_04, pOpt_10, pOpt_25, pOpt_25NW, pOpt_40, pOpt_40NW ])
pErrs = n.array([pCov_04.diagonal()**0.5, pCov_10.diagonal()**0.5, pCov_25.diagonal()**0.5, pCov_25NW.diagonal()**0.5, pCov_40.diagonal()**0.5, pCov_40NW.diagonal()**0.5 ])
print "best params=", pOpts.T
print "err=", pErrs.T
