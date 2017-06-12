import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
import lib_functions_1pt as lib
import os
import sys

fw = open(join(os.environ['MVIR_DIR'],"HMF_fit.txt"), 'w')

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
from scipy.stats import chi2 as stc2

# plotting modules
import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=14
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
dn_L04 = lambda sigma, a, p, q :  (((bias(sigma, a, p, q) * lib.hmf.growth_factor)**2. * (lib.covariance_factor[0]))**2. + lib.shot_simple(sigma, 400.**3.)  )**0.5
dn_L10 = lambda sigma, a, p, q :  (((bias(sigma, a, p, q) * lib.hmf.growth_factor)**2. * (lib.covariance_factor[1]) )**2. + lib.shot_simple(sigma, 1000.**3.) )**0.5  
dn_L25 = lambda sigma, a, p, q : (((bias(sigma, a, p, q) * lib.hmf.growth_factor)**2. * (lib.covariance_factor[2]) )**2. + lib.shot_simple(sigma, 2500.**3.) )**0.5  
dn_L40 = lambda sigma, a, p, q :  (((bias(sigma, a, p, q) * lib.hmf.growth_factor)**2. * (lib.covariance_factor[3]) )**2. + lib.shot_simple(sigma, 4000.**3.) )**0.5  
dn_L80 = lambda sigma, a, p, q :  (((bias(sigma, a, p, q) * lib.hmf.growth_factor)**2. * (lib.covariance_factor[3]) )**2. + lib.shot_simple(sigma, 8000.**3.) )**0.5  

# off diagonal error 
dn_cov_L04 = lambda s1, s2, a, p, q : (dn_L04(s1, a, p, q)*dn_L04(s2, a, p, q))**0.5 
dn_cov_L10 = lambda s1, s2, a, p, q : (dn_L10(s1, a, p, q)*dn_L10(s2, a, p, q) )**0.5
dn_cov_L25 = lambda s1, s2, a, p, q : (dn_L25(s1, a, p, q)*dn_L25(s2, a, p, q) )**0.5
dn_cov_L40 = lambda s1, s2, a, p, q : (dn_L40(s1, a, p, q)*dn_L40(s2, a, p, q) )**0.5
dn_cov_L80 = lambda s1, s2, a, p, q : (dn_L80(s1, a, p, q)*dn_L80(s2, a, p, q) )**0.5

# opens the data 

#Quantity studied
qty = "mvir"
version = 'v4'
fw.write(version + '\n')
# working directory
dir = join(os.environ['MVIR_DIR'])
# loads summary file
data = fits.open( join(dir, qty+"_summary.fits"))[1].data

NminCount = 1000
logNpmin = 3
fw.write("NminCount="+str(NminCount) + '\n')
fw.write("logNpmin="+str(logNpmin) + '\n')

zmin = -0.01
zmax = 0.001
fw.write("zmin="+str(zmin) + '\n')
fw.write("zmax="+str(zmax) + '\n')

# x coordinates definition
logsig = -n.log10(data['sigmaM'])#
lognu = n.log10(data['nu2']**0.5)
#log_mvir = data["log_"+qty]
log_mvir = data["log_"+qty] - n.log10(cosmo.h)
mvir = 10**data["log_"+qty] / cosmo.h



print "------------------------------------------"
print "------------------------------------------"
print "------------------------------------------"
print "------            1d central              -----"
print "------------------------------------------"
print "------------------------------------------"
print "------------------------------------------"

#=======================
#=======================
cos = 'cen'
fw.write("c o s="+cos+ '\n')

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
ok1 = (zSel) & (mSel) & (mSel2) & (nSelCen)
# selection per box :
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')
DS80=(data["boxName"]=='DS_8Gpc')

print "ST01 fit MULTIDARK ----------------------"
fw.write('\n')
fw.write('ST01 fit MULTIDARK \n')

ok = (ok1) & (DS80==False)
x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_cen"][ok]**2. + data["dN_counts_cen"][ok]**(-1.))**(0.5)

ps = n.array([0.333, 0.794, 0.247])
log_fsigma = lambda logsigma, A, a, p : n.log10(lib.f_BH(10**-logsigma, A, a, p, 1.))

pOpt, pCov=curve_fit(log_fsigma, x_data, y_data, ps, sigma = y_data_err, maxfev=50000000)#, bounds=boundaries)
chi2 = n.sum(((log_fsigma(x_data, pOpt[0], pOpt[1], pOpt[2])-y_data)/y_data_err)**2. ) 
ndof = (len(x_data) - len(ps)) 
print "best params=", pOpt
print "err=", pCov.diagonal()**0.5
print "chi2 ", chi2, ndof, chi2/ndof
print "P chi2 1-cdf", 1-stc2.cdf(int(chi2),ndof)
print "---------------------------------------------------"
fw.write("best params = "+ str(n.round(pOpt,5)) + '\n')
fw.write("err params = "+ str(n.round(pCov.diagonal()**0.5,5)) + '\n')
fw.write("chi2 params = "+ str(n.round(chi2,2))+", "+str( ndof)+", "+str(n.round(chi2/ndof,2)) + '\n')
fw.write("P chi2 1-cdf = "+ str(n.round(1-stc2.cdf(int(chi2),ndof),5)) + '\n')

pOpt_ST01 = pOpt
pErr_ST01 = pCov.diagonal()**0.5

n.savetxt(join(os.environ['MVIR_DIR'],"mvirFunction_parameters_ST01_MD_MFonly_fit.txt"), n.transpose([pOpt_ST01, pErr_ST01]), header="A a p", fmt='%s')

print "ST01 fit DARKSKIES ----------------------"
fw.write('\n')
fw.write('ST01 fit DARKSKIES \n')

ok = (ok1) & (DS80)

x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_cen"][ok]**2. + data["dN_counts_cen"][ok]**(-1.))**(0.5)

ps = n.array([0.333, 0.794, 0.247])
log_fsigma = lambda logsigma, A, a, p : n.log10(lib.f_BH(10**-logsigma, A, a, p, 1.))

pOpt, pCov=curve_fit(log_fsigma, x_data, y_data, ps, sigma = y_data_err, maxfev=50000000)#, bounds=boundaries)
chi2 = n.sum(((log_fsigma(x_data, pOpt[0], pOpt[1], pOpt[2])-y_data)/y_data_err)**2. ) 
ndof = (len(x_data) - len(ps)) 
print "best params=", pOpt
print "err=", pCov.diagonal()**0.5
print "chi2 ", chi2, ndof, chi2/ndof
print "P chi2 1-cdf", 1-stc2.cdf(int(chi2),ndof)
print "---------------------------------------------------"
fw.write("best params = "+ str(n.round(pOpt,5)) + '\n')
fw.write("err params = "+ str(n.round(pCov.diagonal()**0.5,5)) + '\n')
fw.write("chi2 params = "+ str(n.round(chi2,2))+", "+str( ndof)+", "+str(n.round(chi2/ndof,2)) + '\n')
fw.write("P chi2 1-cdf = "+ str(n.round(1-stc2.cdf(int(chi2),ndof),5)) + '\n')
pOpt_ST01 = pOpt
pErr_ST01 = pCov.diagonal()**0.5

n.savetxt(join(os.environ['MVIR_DIR'],"mvirFunction_parameters_ST01_DS_MFonly_fit.txt"), n.transpose([pOpt_ST01, pErr_ST01]), header="A a p", fmt='%s')

print "---------------------------------------------------"
print "---------------------------------------------------"
print "---------------------------------------------------"

print "BATT 2011 MULTIDARK"
fw.write('\n')
fw.write('BATT fit MULTIDARK \n')

ok = (ok1) & (DS80==False)

x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_cen"][ok]**2. + data["dN_counts_cen"][ok]**(-1.))**(0.5)

ps = n.array([A0, a0, p0, q0])
log_fsigma = lambda logsigma, A, a, p, q : n.log10(lib.f_BH(10**-logsigma, A, a, p, q))

pOpt, pCov=curve_fit(log_fsigma, x_data, y_data, ps, sigma = y_data_err, maxfev=50000000)#, bounds=boundaries)
chi2 = n.sum(((log_fsigma(x_data, pOpt[0], pOpt[1], pOpt[2], pOpt[3])-y_data)/y_data_err)**2. ) 
ndof = (len(x_data) - len(ps)) 
print "best params=", pOpt
print "err=", pCov.diagonal()**0.5
print "chi2 ", chi2, ndof, chi2/ndof
print "P chi2 1-cdf", 1-stc2.cdf(int(chi2),ndof)
print "---------------------------------------------------"
fw.write("best params = "+ str(n.round(pOpt,5)) + '\n')
fw.write("err params = "+ str(n.round(pCov.diagonal()**0.5,5)) + '\n')
fw.write("chi2 params = "+ str(n.round(chi2,2))+", "+str( ndof)+", "+str(n.round(chi2/ndof,2)) + '\n')
fw.write("P chi2 1-cdf = "+ str(n.round(1-stc2.cdf(int(chi2),ndof),5)) + '\n')
A1, a1, p1, q1 = pOpt
A1_err, a1_err, p1_err, q1_err = pCov.diagonal()**0.5

n.savetxt(join(os.environ['MVIR_DIR'],"mvirFunction_parameters_BA11_MD_MFonly_fit.txt"), n.transpose([pOpt, pCov.diagonal()**0.5]), header="A a p q", fmt='%s')

print "BATT 2011 DARKSKIES"
fw.write('\n')
fw.write('BATT fit DARKSKIES \n')

ok = (ok1) & (DS80)

x_data = logsig[ok]
y_data = log_MF[ok]
y_data_err = (data["std90_pc_cen"][ok]**2. + data["dN_counts_cen"][ok]**(-1.))**(0.5)

ps = n.array([A0, a0, p0, q0])
log_fsigma = lambda logsigma, A, a, p, q : n.log10(lib.f_BH(10**-logsigma, A, a, p, q))

pOpt, pCov=curve_fit(log_fsigma, x_data, y_data, ps, sigma = y_data_err, maxfev=50000000)#, bounds=boundaries)
chi2 = n.sum(((log_fsigma(x_data, pOpt[0], pOpt[1], pOpt[2], pOpt[3])-y_data)/y_data_err)**2. ) 
ndof = (len(x_data) - len(ps)) 
print "best params=", pOpt
print "err=", pCov.diagonal()**0.5
print "chi2 ", chi2, ndof, chi2/ndof
print "P chi2 1-cdf", 1-stc2.cdf(int(chi2),ndof)
print "---------------------------------------------------"
fw.write("best params = "+ str(n.round(pOpt,5)) + '\n')
fw.write("err params = "+ str(n.round(pCov.diagonal()**0.5,5)) + '\n')
fw.write("chi2 params = "+ str(n.round(chi2,2))+", "+str( ndof)+", "+str(n.round(chi2/ndof,2)) + '\n')
fw.write("P chi2 1-cdf = "+ str(n.round(1-stc2.cdf(int(chi2),ndof),5)) + '\n')
A1, a1, p1, q1 = pOpt
A1_err, a1_err, p1_err, q1_err = pCov.diagonal()**0.5

n.savetxt(join(os.environ['MVIR_DIR'],"mvirFunction_parameters_BA11_DS_MFonly_fit.txt"), n.transpose([pOpt, pCov.diagonal()**0.5]), header="A a p q", fmt='%s')

ok = (ok1) #& (DS80==False)

def plotSel(MDsel, label):
	x_data = logsig[ok & MDsel]
	y_data = log_MF[ok & MDsel]
	y_data_err = (data["std90_pc_cen"][ok & MDsel]**2. + data["dN_counts_cen"][ok & MDsel]**(-1.))**(0.5)
	y_model = log_fsigma(x_data, pOpt[0], pOpt[1], pOpt[2], pOpt[3])
	if len(x_data)>0:
		p.errorbar(x_data, 10**(y_data-y_model), yerr = y_data_err , rasterized=True, fmt='none', label=label)


p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
plotSel(MD04, "M04")
plotSel(MD10, "M10")
plotSel(MD25, "M25")
plotSel(MD40, "M40")
plotSel(MD25NW, "M25n")
plotSel(MD40NW, "M40n")
plotSel(DS80, "D80")
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

		
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
plotAll(MD04, "M04")
plotAll(MD10, "M10")
plotAll(MD25, "M25")
plotAll(MD25NW, "M25n")
plotAll(MD40, "M40")
plotAll(MD40NW, "M40n")
plotAll(DS80, "D80")

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

fw.close()