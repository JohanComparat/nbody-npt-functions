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
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p

# mass function theory
from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

#lib.covariance_factor
#lib.f_BH(sigma, 0.333, 0.788, 0.807, 1.795)
bias = lambda sigma : lib.b_BH(sigma, a=0.8915, p=0.5524, q=1.578)

# sample variance equation 16
var_sv_L04 = lambda sa, sb, fc=lib.covariance_factor_jk : bias(sa) * bias(sb) * (lib.hmf.growth_factor)**2. * fc[0]
var_sv_L10 = lambda sa, sb, fc=lib.covariance_factor_jk : bias(sa) * bias(sb) * (lib.hmf.growth_factor)**2. * fc[1]
var_sv_L25 = lambda sa, sb, fc=lib.covariance_factor_jk : bias(sa) * bias(sb) * (lib.hmf.growth_factor)**2. * fc[2]
var_sv_L40 = lambda sa, sb, fc=lib.covariance_factor_jk : bias(sa) * bias(sb) * (lib.hmf.growth_factor)**2. * fc[3]
var_sv_L80 = lambda sa, sb, fc=lib.covariance_factor_jk : bias(sa) * bias(sb) * (lib.hmf.growth_factor)**2. * fc[4]

# shot noise variance 
var_sn_L04 = lambda sigma, lbox=40.   : lib.shot_simple(sigma, lbox**3.)
var_sn_L10 = lambda sigma, lbox=100. : lib.shot_simple(sigma, lbox**3.)
var_sn_L25 = lambda sigma, lbox=250. : lib.shot_simple(sigma, lbox**3.)
var_sn_L40 = lambda sigma, lbox=400. : lib.shot_simple(sigma, lbox**3.)
var_sn_L80 = lambda sigma, lbox=800. : lib.shot_simple(sigma, lbox**3.)

# total covariance equation 17

print '--------------------------------------------------'
print '--------------------------------------------------'
print '------------------MATRIX------------------'
print '--------------------------------------------------'
print '--------------------------------------------------'

#Quantity studied
qty = "mvir"
version = 'v3'
# measurement files
fileC = n.array(glob.glob( join(os.environ['MD_DIR'], "MD_*Gpc*",  version, qty,"out_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MD_DIR'], "MD_*Gpc*", version, qty,"out_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MD_DIR'], "MD_*Gpc*", version, qty,"out_*_Satellite_JKresampling.pkl")))

# redshift 0 data
iis = [8, 13, 29, 31, 60, -10]#[-1, -2, -4, -9, -22, 3]
print fileC[iis]

def plot_COV(fileC, binFile, gt14):
	"""
	Reads the data and extract mass function covariance matrix
	
	First gets the basic info about the data
	
	Then creates a dimension-less mass function matrix (1000 mass functions)
	
	"""
	# hig mass end
	if gt14:
		snFactor = n.array([4.8, 5.5, 8., 11., 12.])
		svFactor = n.array([1./2., 1./3.2, 4.5, 5, 4.5])
		binW = 0.025 
		plotDir = "covariance_gt14"
	#low mass end
	else:
		snFactor = n.array([4.8, 5.5, 8., 11., 12])
		svFactor = n.array([1./2., 1./3.2, 4.5, 5, 4.5])
		binW = 0.125 
		plotDir = "covariance_lt14"

	ddd, boxName = lib.convert_pkl_massFunction_covarianceMatrix(fileC, binFile, qty='mvir', delta_wrt='mean', gt14=gt14)
	#print ddd
	print len(ddd)
	f_mean, f_matrix, count_matrix, sigma, mass = ddd
	cv = n.cov(f_matrix.T, ddof=0)#/(dxcv*dycv)
	ctotal = (cv + count_matrix**(-0.5) )#/(dxcv*dycv)

	var_total_sv_L04 = lambda s1, s2 : var_sv_L04(s1, s2, lib.covariance_factor_jk2) / svFactor[0]
	var_total_sv_L10 = lambda s1, s2 : var_sv_L10(s1, s2, lib.covariance_factor_jk2) / svFactor[1]
	var_total_sv_L25 = lambda s1, s2 : var_sv_L25(s1, s2, lib.covariance_factor_jk2) / svFactor[2]
	var_total_sv_L40 = lambda s1, s2 : var_sv_L40(s1, s2, lib.covariance_factor_jk2) / svFactor[3]
	var_total_sv_L80 = lambda s1, s2 : var_sv_L80(s1, s2, lib.covariance_factor_jk2) / svFactor[4]

	var_total_sn_L04 = lambda s1, s2 : lib.shot_double(s1, s2, (  40.)**3.  , binW) / snFactor[0] 
	var_total_sn_L10 = lambda s1, s2 : lib.shot_double(s1, s2, ( 100.)**3. , binW) / snFactor[1] 
	var_total_sn_L25 = lambda s1, s2 : lib.shot_double(s1, s2, ( 250.)**3. , binW) / snFactor[2] 
	var_total_sn_L40 = lambda s1, s2 : lib.shot_double(s1, s2, ( 400.)**3. , binW) / snFactor[3] 
	var_total_sn_L80 = lambda s1, s2 : lib.shot_double(s1, s2, ( 800.)**3. , binW) / snFactor[4] 

	var_total_L04 = lambda s1, s2 : (var_total_sn_L04(s1, s2) + var_total_sv_L04(s1, s2))
	var_total_L10 = lambda s1, s2 : (var_total_sn_L10(s1, s2) + var_total_sv_L10(s1, s2))
	var_total_L25 = lambda s1, s2 : (var_total_sn_L25(s1, s2) + var_total_sv_L25(s1, s2))
	var_total_L40 = lambda s1, s2 : (var_total_sn_L40(s1, s2) + var_total_sv_L40(s1, s2))
	var_total_L80 = lambda s1, s2 : (var_total_sn_L80(s1, s2) + var_total_sv_L80(s1, s2))
	
	xcv, ycv = n.meshgrid(sigma, sigma)
	
	model = {"MD_0.4Gpc": var_total_L04, "MD_1Gpc": var_total_L10, "MD_2.5Gpc": var_total_L25, "MD_4Gpc": var_total_L40,"MD_2.5GpcNW": var_total_L25, "MD_4GpcNW": var_total_L40 , "DS_8Gpc": var_total_L80 }
	model_sn = {"MD_0.4Gpc": var_total_sn_L04, "MD_1Gpc": var_total_sn_L10, "MD_2.5Gpc": var_total_sn_L25, "MD_4Gpc": var_total_sn_L40,"MD_2.5GpcNW": var_total_sn_L25, "MD_4GpcNW": var_total_sn_L40 , "DS_8Gpc": var_total_sn_L80 }
	model_sv = {"MD_0.4Gpc": var_total_sv_L04, "MD_1Gpc": var_total_sv_L10, "MD_2.5Gpc": var_total_sv_L25, "MD_4Gpc": var_total_sv_L40,"MD_2.5GpcNW": var_total_sv_L25, "MD_4GpcNW": var_total_sv_L40 , "DS_8Gpc": var_total_sv_L80 }

	nickname = {"MD_0.4Gpc": "M04", "MD_1Gpc": "M10", "MD_2.5Gpc": "M25", "MD_4Gpc": "M40","MD_2.5GpcNW": "M25n", "MD_4GpcNW": "M40n" , "DS_8Gpc": "D80" }
	p.figure(1,(6,6))
	p.axes([0.17, 0.17, 0.78, 0.78])
	p.title(nickname[boxName])#+" total")
	tp=(xcv>ycv)# (xcv!=ycv)
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=n.log10(ctotal[tp]), s=25, edgecolors='none',vmin=-2,vmax=2., marker='s', label='data')
	tp=(xcv<ycv)# (xcv!=ycv)
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=n.log10(model[boxName](xcv, ycv)[tp]), s=25, edgecolors='none', vmin=-2,vmax=2, marker='o', label='model')
	cb = p.colorbar(shrink=0.7)
	cb.set_label(r'$\log_{10}(C(\sigma_1, \sigma_2))$')
	p.xlabel(r'$log_{10}(\sigma_1^{-1})$')
	p.ylabel(r'$log_{10}(\sigma_2^{-1})$')
	#p.xlim((-0.7, 0.6))
	#p.ylim((-0.7, 0.6))
	p.grid()
	gl=p.legend(loc=0, frameon=False)
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_"+boxName+"_dataTOT.png"))
	p.clf()
	"""
	p.figure(1,(6,6))
	p.axes([0.17, 0.17, 0.78, 0.78])
	p.title(nickname[boxName]+" sn")
	tp=(xcv>ycv)# (xcv!=ycv)
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=n.log10(ctotal[tp]), s=15, edgecolors='none',vmin=-2,vmax=2., marker='s')
	tp=(xcv<ycv)# (xcv!=ycv)
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=n.log10(model_sn[boxName](xcv, ycv)[tp]), s=15, edgecolors='none', vmin=-2,vmax=2, marker='s')
	cb = p.colorbar(shrink=0.7)
	cb.set_label(r'$\log_{10}(C(\sigma_1, \sigma_2))$')
	p.xlabel(r'$log_{10}(\sigma_1^{-1})$')
	p.ylabel(r'$log_{10}(\sigma_2^{-1})$')
	p.xlim((-0.7, 0.6))
	p.ylim((-0.7, 0.6))
	p.grid()
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_"+boxName+"_dataSN.png"))
	p.clf()
	
	p.figure(1,(6,6))
	p.axes([0.17, 0.17, 0.78, 0.78])
	p.title(nickname[boxName]+" sv")
	tp=(xcv>ycv)# (xcv!=ycv)
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=n.log10(ctotal[tp]), s=15, edgecolors='none',vmin=-2,vmax=2., marker='s')
	tp=(xcv<ycv)# (xcv!=ycv)
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=n.log10(model_sv[boxName](xcv, ycv)[tp]), s=15, edgecolors='none', vmin=-2,vmax=2, marker='s')
	cb = p.colorbar(shrink=0.7)
	cb.set_label(r'$\log_{10}(C(\sigma_1, \sigma_2))$')
	p.xlabel(r'$log_{10}(\sigma_1^{-1})$')
	p.ylabel(r'$log_{10}(\sigma_2^{-1})$')
	p.xlim((-0.7, 0.6))
	p.ylim((-0.7, 0.6))
	p.grid()
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_"+boxName+"_dataSV.png"))
	p.clf()
	"""
	p.figure(2,(6,6))
	p.axes([0.17, 0.17, 0.78, 0.78])
	p.title(nickname[boxName])
	tp=(xcv>ycv)
	residuals = n.log10(ctotal[tp]/model[boxName](xcv, ycv)[tp])
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=residuals, s=25, edgecolors='none', vmin=-0.5,vmax=0.5, marker='s')
	cb = p.colorbar(shrink=0.7)
	cb.set_label(r'$\log_{10}$(data/model)')
	#p.xlim((-0.7, 0.6))
	p.grid()
	#p.ylim((-0.7, 0.6))
	p.xlabel(r'$log_{10}(\sigma_1^{-1})$')
	p.ylabel(r'$log_{10}(\sigma_2^{-1})$')
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_"+boxName+"_model.png"))
	p.clf()
	
	p.figure(12)
	for el in f_matrix:
		p.plot(-n.log10(sigma), n.log10(el), 'k', alpha=0.1)
		
	p.plot(-n.log10(sigma), n.log10(f_mean),'r',lw=2)
	p.ylim((-3., -0.4))
	p.xlim((-0.5, 0.5))
	p.xlabel(r'$log_{10}(\sigma^{-1})$')
	p.ylabel(r'$\log_{10}\left[ \frac{M}{\rho_m} \frac{dn}{d\ln M} \left|\frac{d\ln M }{d\ln \sigma}\right|\right] $') 
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"cv_test"+boxName+".png"))
	p.clf()
	
	p.figure(2,(6,6))
	p.axes([0.17, 0.17, 0.78, 0.78])
	p.title(nickname[boxName])
	p.errorbar(-n.log10(xcv.diagonal()), n.log10(ctotal.diagonal()),yerr=0.1, label='data')
	y=model[boxName](xcv, ycv)
	p.plot(-n.log10(xcv.diagonal()), n.log10(y.diagonal()), 'k', label='model')
	y=model_sn[boxName](xcv, ycv)
	p.plot(-n.log10(xcv.diagonal()), n.log10(y.diagonal()), 'k-.', label='model sn')
	y=model_sv[boxName](xcv, ycv)
	p.plot(-n.log10(xcv.diagonal()), n.log10(y.diagonal()), 'k--', label='model sv')
	#p.xlim((-0.7, 0.6))
	p.grid()
	#p.ylim((-0.7, 0.6))
	p.xlabel(r'$log_{10}(\sigma^{-1})$')
	p.ylabel(r'$log_{10}(C_{Diag})$')
	gl=p.legend(loc=0, frameon=False)
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_"+boxName+"_diagonal.png"))
	p.clf()
	
	p.figure(2,(6,6))
	p.axes([0.17, 0.17, 0.78, 0.78])
	p.title(nickname[boxName])
	y=model[boxName](xcv, ycv).diagonal() / ctotal.diagonal()
	p.plot(-n.log10(xcv.diagonal()), y, 'k', label='model')
	#p.xlim((-0.7, 0.6))
	p.grid()
	#p.ylim((-0.7, 0.6))
	p.xlabel(r'$log_{10}(\sigma^{-1})$')
	p.ylabel(r'diagonal $C_{model}/C_{data}$')
	gl=p.legend(loc=0, frameon=False)
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_"+boxName+"_diagonal_residual.png"))
	p.clf()
	
	
	return f_mean, f_matrix, count_matrix, sigma, mass, residuals

def compCov(gt):
	dout = []
	if gt:
		plotDir = "covariance_gt14"
		for ii in iis:
			dout.append( plot_COV(fileC[ii], fileB[ii], gt ))
	#low mass end
	else:
		plotDir = "covariance_lt14"
		for ii in iis[:-2]:
			dout.append( plot_COV(fileC[ii], fileB[ii], gt ))

	"""
	fileC = n.array(glob.glob( join(os.environ['DS_DIR'], version, qty,"ds*_Central_JKresampling.pkl")))
	fileB = n.array(glob.glob( join( os.environ['DS_DIR'], version, qty,"ds*_"+qty+"_JKresampling.bins")))
	fileS = n.array(glob.glob( join( os.environ['DS_DIR'], version, qty,"ds*_Satellite_JKresampling.pkl")))
	plot_COV(fileC[0], fileB[0])
	ii=0
	dout.append( plot_COV(fileC[ii], fileB[ii]) )
	"""
	f_mean, f_matrix, count_matrix, sigma, mass, residuals = n.transpose(dout)

	MDnames= n.array(['M04', 'M10', 'M25','M25n','M40','M40n', 'D80'])

	xmax = n.empty(len(MDnames))
	p.figure(6,(6,6))
	for ii,el in enumerate(residuals):
		triM = n.tril(el, k=-1)
		nn,bb,pp=p.hist(triM[(triM>-2000)&(triM<2000)&(triM!=0)], bins=25, histtype='step',label=MDnames[ii], normed=True)
		xmax[ii]=(bb[n.argmax(nn)]+bb[n.argmax(nn)+1])/2.
		
	p.xlabel('log(data/model)')
	p.ylabel('normed counts')
	p.title('Off diagonal covariance')
	p.grid()
	p.xlim((-0.4, 0.4))
	gl=p.legend(loc=0, frameon=False)
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"residuals.png"))
	p.clf()

	print "data/model most populated bin"
	for nm, xm in  zip(MDnames, xmax)[:-1]:
		print nm, 10**xm

compCov(True)
compCov(False)

sys.exit()


print '--------------------------------------------------'
print '--------------------------------------------------'
print '------------------DIAGONAL------------------'
print '--------------------------------------------------'
print '--------------------------------------------------'

# projection for the diagonal error using 90% of the volume: only one shot noise term
var_sv_L04_90 = lambda sa, sb : var_sv_L04(sa, sb, lib.covariance_factor_90)
var_sv_L10_90 = lambda sa, sb : var_sv_L10(sa, sb, lib.covariance_factor_90)
var_sv_L25_90 = lambda sa, sb : var_sv_L25(sa, sb, lib.covariance_factor_90)
var_sv_L40_90 = lambda sa, sb : var_sv_L40(sa, sb, lib.covariance_factor_90)
var_sv_L80_90 = lambda sa, sb : var_sv_L80(sa, sb, lib.covariance_factor_90)

var_sn_L04_90 = lambda sigma:  var_sn_L04(sigma, lbox=400.*0.9**(1./3.)   )
var_sn_L10_90 = lambda sigma : var_sn_L10(sigma, lbox=1000.*0.9**(1./3.) )
var_sn_L25_90 = lambda sigma : var_sn_L25(sigma, lbox=2500.*0.9**(1./3.) ) 
var_sn_L40_90 = lambda sigma : var_sn_L40(sigma, lbox=4000.*0.9**(1./3.) ) 
var_sn_L80_90 = lambda sigma : var_sn_L80(sigma, lbox=8000.*0.9**(1./3.) ) 

var_total_L04_90 = lambda s1, s2 : var_sn_L04_90(s1) + var_sv_L04_90(s1, s2) 
var_total_L10_90 = lambda s1, s2 : var_sn_L10_90(s1) + var_sv_L10_90(s1, s2) 
var_total_L25_90 = lambda s1, s2 : var_sn_L25_90(s1) + var_sv_L25_90(s1, s2) 
var_total_L40_90 = lambda s1, s2 : var_sn_L40_90(s1) + var_sv_L40_90(s1, s2) 
var_total_L80_90 = lambda s1, s2 : var_sn_L80_90(s1) + var_sv_L80_90(s1, s2) 


qty = 'mvir'
dir = join(os.environ['MVIR_DIR'])
# loads summary file
dataMF = fits.open( join(dir, qty+"_summary.fits"))[1].data
zzero = (dataMF['redshift']==0) & (dataMF['log_mvir']>3+dataMF['logMpart']) & (dataMF['dN_counts_cen'] > 10 )
dlnSigM = abs(n.log(dataMF['log_mvir_max']-dataMF['log_mvir_min'])*dataMF['dlnsigmaMdlnM'])

p.figure(0, (6,6))
p.axes([0.17, 0.17, 0.78, 0.78])
x = n.logspace(-0.51, 0.51, 25)

zSel = (zzero) & (dataMF['boxName'] =='MD_0.4Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["std90_pc_cen"][zSel], 'b+')
p.plot(-n.log10(x), var_sv_L04_90(x,x)**0.5, 'b--')
p.plot(-n.log10(x), var_sn_L04_90(x)**0.5, 'b', ls='dotted')
p.plot(-n.log10(x), var_total_L04_90(x,x)**0.5, 'b')

zSel = (zzero) & (dataMF['boxName'] =='MD_1Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["std90_pc_cen"][zSel], 'g+')
p.plot(-n.log10(x), var_sv_L10_90(x,x)**0.5, 'g--')
p.plot(-n.log10(x), var_sn_L10_90(x)**0.5, 'g', ls='dotted')
p.plot(-n.log10(x), var_total_L10_90(x,x)**0.5, 'g')


zSel = (zzero) & (dataMF['boxName'] =='MD_2.5Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["std90_pc_cen"][zSel], 'r+')
p.plot(-n.log10(x), var_sv_L25_90(x,x)**0.5, 'r--')
p.plot(-n.log10(x), var_sn_L25_90(x)**0.5, 'r', ls='dotted')
p.plot(-n.log10(x), var_total_L25_90(x,x)**0.5, 'r')

zSel = (zzero) & (dataMF['boxName'] =='MD_4Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["std90_pc_cen"][zSel], 'k+')
p.plot(-n.log10(x), var_sv_L40_90(x,x)**0.5, 'k--')
p.plot(-n.log10(x), var_sn_L40_90(x)**0.5, 'k', ls='dotted')
p.plot(-n.log10(x), var_total_L40_90(x,x)**0.5, 'k')

zSel = (zzero) & (dataMF['boxName'] =='DS_8Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["std90_pc_cen"][zSel], 'm+')
p.plot(-n.log10(x), var_sv_L80_90(x,x)**0.5, 'm--')
p.plot(-n.log10(x), var_sn_L80_90(x)**0.5, 'm', ls='dotted')
p.plot(-n.log10(x), var_total_L80_90(x,x)**0.5, 'm')

p.plot(0.,0., 'k', label='sample variance', ls='dashed')
p.plot(0.,0., 'k', label='shot noise', ls='dotted')
p.plot(0.,0., 'k', label='sum', ls='solid')
# p.xscale('log')
p.yscale('log')
p.ylim((2e-4, 10))
p.xlim((-0.5, 0.7))
p.grid()
p.xlabel(r'$log_{10}(\sigma^{-1})$')
p.ylabel(r'fractional error $\sqrt{C_{model}(\sigma,\sigma)}$')
gl=p.legend(loc=0, frameon=False)
p.title('jackknife')
#gl.set_frame_on(False)
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"dn-sigma-jackknife.png"))
p.clf()


p.figure(0, (6,6))
p.axes([0.17, 0.17, 0.78, 0.78])
x = n.logspace(-0.51, 0.51, 25)

zSel = (zzero) & (dataMF['boxName'] =='MD_0.4Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["dN_counts_cen"][zSel]**(-0.5), 'bx', label='M04')
p.plot(-n.log10(x), var_sv_L04_90(x,x)**0.5, 'b--')
p.plot(-n.log10(x), var_sn_L04_90(x)**0.5, 'b', ls='dotted')
p.plot(-n.log10(x), var_total_L04_90(x,x)**0.5, 'b')

zSel = (zzero) & (dataMF['boxName'] =='MD_1Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["dN_counts_cen"][zSel]**(-0.5), 'gx', label='M10')
p.plot(-n.log10(x), var_sv_L10_90(x,x)**0.5, 'g--')
p.plot(-n.log10(x), var_sn_L10_90(x)**0.5, 'g', ls='dotted')
p.plot(-n.log10(x), var_total_L10_90(x,x)**0.5, 'g')

zSel = (zzero) & (dataMF['boxName'] =='MD_2.5Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["dN_counts_cen"][zSel]**(-0.5), 'rx', label='M25')
p.plot(-n.log10(x), var_sv_L25_90(x,x)**0.5, 'r--')
p.plot(-n.log10(x), var_sn_L25_90(x)**0.5, 'r', ls='dotted')
p.plot(-n.log10(x), var_total_L25_90(x,x)**0.5, 'r')

zSel = (zzero) & (dataMF['boxName'] =='MD_4Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["dN_counts_cen"][zSel]**(-0.5), 'kx', label='M40')
p.plot(-n.log10(x), var_sv_L40_90(x,x)**0.5, 'k--')
p.plot(-n.log10(x), var_sn_L40_90(x)**0.5, 'k', ls='dotted')
p.plot(-n.log10(x), var_total_L40_90(x,x)**0.5, 'k')

zSel = (zzero) & (dataMF['boxName'] =='DS_8Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), dataMF["dN_counts_cen"][zSel]**(-0.5), 'mx', label='D80')
p.plot(-n.log10(x), var_sv_L80_90(x,x)**0.5, 'm--')
p.plot(-n.log10(x), var_sn_L80_90(x)**0.5, 'm', ls='dotted')
p.plot(-n.log10(x), var_total_L80_90(x,x)**0.5, 'm')

p.yscale('log')
p.ylim((2e-4, 10))
p.xlim((-0.5, 0.7))
p.grid()
p.xlabel(r'$log_{10}(\sigma^{-1})$')
p.ylabel(r'fractional error $\sqrt{C_{model}(\sigma,\sigma)}$')
gl=p.legend(loc=0, frameon=False)
p.title('poisson')
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"dn-sigma-poisson.png"))
p.clf()

p.figure(0, (6,6))
p.axes([0.17, 0.17, 0.78, 0.78])
x = n.logspace(-0.51, 0.51, 25)

zSel = (zzero) & (dataMF['boxName'] =='MD_0.4Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), (dataMF["std90_pc_cen"][zSel]**2. + dataMF["dN_counts_cen"][zSel]**(-1.))**(0.5), 'bx', label='MD04')
p.plot(-n.log10(x), var_sv_L04_90(x,x)**0.5, 'b--')
p.plot(-n.log10(x), var_sn_L04_90(x)**0.5, 'b', ls='dotted')
p.plot(-n.log10(x), var_total_L04_90(x,x)**0.5, 'b')

zSel = (zzero) & (dataMF['boxName'] =='MD_1Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), (dataMF["std90_pc_cen"][zSel]**2. + dataMF["dN_counts_cen"][zSel]**(-1.))**(0.5), 'gx', label='MD10')
p.plot(-n.log10(x), var_sv_L10_90(x,x)**0.5, 'g--')
p.plot(-n.log10(x), var_sn_L10_90(x)**0.5, 'g', ls='dotted')
p.plot(-n.log10(x), var_total_L10_90(x,x)**0.5, 'g')

zSel = (zzero) & (dataMF['boxName'] =='MD_2.5Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), (dataMF["std90_pc_cen"][zSel]**2. + dataMF["dN_counts_cen"][zSel]**(-1.))**(0.5), 'rx', label='MD25')
p.plot(-n.log10(x), var_sv_L25_90(x,x)**0.5, 'r--')
p.plot(-n.log10(x), var_sn_L25_90(x)**0.5, 'r', ls='dotted')
p.plot(-n.log10(x), var_total_L25_90(x,x)**0.5, 'r')

zSel = (zzero) & (dataMF['boxName'] =='MD_4Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), (dataMF["std90_pc_cen"][zSel]**2. + dataMF["dN_counts_cen"][zSel]**(-1.))**(0.5), 'kx', label='MD40')
p.plot(-n.log10(x), var_sv_L40_90(x,x)**0.5, 'k--')
p.plot(-n.log10(x), var_sn_L40_90(x)**0.5, 'k', ls='dotted')
p.plot(-n.log10(x), var_total_L40_90(x,x)**0.5, 'k')

zSel = (zzero) & (dataMF['boxName'] =='DS_8Gpc')
p.plot(-n.log10(dataMF['sigmaM'][zSel]), (dataMF["std90_pc_cen"][zSel]**2. + dataMF["dN_counts_cen"][zSel]**(-1.))**(0.5), 'mx', label='DS80')
p.plot(-n.log10(x), var_sv_L80_90(x,x)**0.5, 'm--')
p.plot(-n.log10(x), var_sn_L80_90(x)**0.5, 'm', ls='dotted')
p.plot(-n.log10(x), var_total_L80_90(x,x)**0.5, 'm')

#p.plot(x, lib.shot_noise(x, 400.**3.)**0.5, label='04sn')
#p.plot(x, lib.shot_noise(x, 1000.**3.)**0.5, label='10sn')
#p.plot(x, lib.shot_noise(x, 2500.**3.)**0.5, label='25sn')
#p.plot(x, lib.shot_noise(x, 4000.**3.)**0.5, label='40sn')
# p.xscale('log')
p.yscale('log')
p.ylim((2e-4, 10))
p.xlim((-0.5, 0.7))
p.grid()
p.xlabel(r'$log_{10}(\sigma^{-1})$')
p.ylabel(r'$\Delta n / n / dln(\sigma)$')
p.ylabel(r'fractional error $\sqrt{C_{model}(\sigma,\sigma)}$')
gl=p.legend(loc=0, frameon=False)
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"dn-sigma-both.png"))
p.clf()

sys.exit()

def plot_CRCoef_mvir(fileC, fileS, binFile):
	boxZN = float(os.path.basename(fileC).split('_')[1])
	print boxZN
	hf, boxLength, boxName, boxRedshift, logmp, boxLengthComoving, massCorrection = lib.get_basic_info(fileC, boxZN, delta_wrt='mean')
	bins = n.log10( 10**n.loadtxt(binFile) * massCorrection )
	logmass = ( bins[1:]  + bins[:-1] )/2.
	mSel = (logmass > 9) & (logmass < 15)
	mass = 10**logmass
	dX = ( 10**bins[1:]  - 10**bins[:-1] )
	dlnbin = (bins[1:]  - bins[:-1])*n.log(10)

	m2sigma = interp1d(hf.M, hf.sigma )
	sigma = m2sigma( mass )
	
	data=cPickle.load(open(fileC,'r'))
	dataS=cPickle.load(open(fileS,'r'))
	
	counts = n.sum(data, axis=0)
	count_matrix = n.outer(counts, counts)/1000.
	cv = (n.cov(data.T, ddof=0)/count_matrix)**0.5
	cr = n.corrcoef(data.T)
	cvS = (n.cov(dataS.T, ddof=0)/count_matrix)**0.5
	crS = n.corrcoef(dataS.T)
	ctotal = cv + count_matrix**(-0.5)
	
	mass2X = interp1d(logmass, n.arange(len(logmass)))
	sigma2X = interp1d(sigma, n.arange(len(logmass)))

	fig = p.figure(0,(6,6))
	mat = p.matshow(cr)
	p.xticks(n.arange(0,len(logmass),5), logmass[n.arange(0,len(logmass),5)],rotation=45)
	p.yticks(n.arange(0,len(logmass),5), logmass[n.arange(0,len(logmass),5)])
	p.axvline(mass2X(logmp+2), lw=2, color='k')
	p.axhline(mass2X(logmp+2), lw=2, color='k')
	#p.axvline(mass2X(logmp+1), lw=2, color='k')
	#p.axhline(mass2X(logmp+1), lw=2, color='k')
	cb = p.colorbar(shrink=0.8)
	cb.set_label(r"R")
	p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
	p.grid()
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"mvir-cr-"+boxName+".png"))
	p.clf()
	
	fig = p.figure(0,(6,6))
	mat = p.matshow(cr)
	p.xticks(n.arange(0,len(sigma),5), n.round(sigma[n.arange(0,len(sigma),5)],3),rotation=45)
	p.yticks(n.arange(0,len(sigma),5), n.round(sigma[n.arange(0,len(sigma),5)],3))
	cb = p.colorbar(shrink=0.8)
	cb.set_label(r"R($\sigma_1, \sigma_2$)")
	p.xlabel(r'$\sigma$')
	p.ylabel(r'$\sigma$')
	p.grid()
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"sigma-cr-"+boxName+".png"))
	p.clf()
	id = int(mass2X(logmp+2.5))
	print id, len(logmass)
	return cr, cv, logmass, sigma, id, ctotal

crs=n.zeros((len(iis),79,79))
cvs=n.zeros((len(iis),79,79))
lgm=n.zeros((len(iis),79))
sigs=n.zeros((len(iis),79))
ctot=n.zeros((len(iis),79,79))
ids=n.zeros((len(iis)))
for ii, el in enumerate(iis):
	crs[ii], cvs[ii], lgm[ii], sigs[ii], ids[ii], ctot[ii] = plot_CRCoef_mvir(fileC[el], fileS[el], fileB[el])


#### SIGMA COVARIANCE

xcr_a, ycr_a, zcr_a = [], [], []
xcv_a, ycv_a, zcv_a = [], [], []
xct_a, yct_a, zct_a = [], [], []

for ii in range(len(ids)):
	xsig_a, ysig_a = n.meshgrid(sigs[ii][ids[ii]:], sigs[ii][ids[ii]:])
	
	zsig_cr_a = crs[ii][ids[ii]:,ids[ii]:]
	zsig_cv_a = cvs[ii][ids[ii]:,ids[ii]:]
	zsig_ct_a = ctot[ii][ids[ii]:,ids[ii]:]

	xsig_b = n.ravel(xsig_a)
	ysig_b = n.ravel(ysig_a)
	zsig_cr_b = n.ravel(zsig_cr_a)
	zsig_cv_b = n.ravel(zsig_cv_a)
	zsig_ct_b = n.ravel(zsig_ct_a)
	
	ok_cr = (n.isnan(zsig_cr_b)==False)&(zsig_cr_b!=n.inf)
	ok_cv = (n.isnan(zsig_cv_b)==False)&(zsig_cv_b!=n.inf)
	ok_ct = (n.isnan(zsig_ct_b)==False)&(zsig_ct_b!=n.inf)
	
	xcr_a.append( xsig_b[ok_cr] )
	ycr_a.append( ysig_b[ok_cr] )
	zcr_a.append( zsig_cr_b[ok_cr] )
	
	xcv_a.append( xsig_b[ok_cv] )
	ycv_a.append( ysig_b[ok_cv] )
	zcv_a.append( zsig_cv_b[ok_cv] )

	xct_a.append( xsig_b[ok_ct] )
	yct_a.append( ysig_b[ok_ct] )
	zct_a.append( zsig_ct_b[ok_ct] )

xcr=n.hstack(xcr_a)
ycr=n.hstack(ycr_a)
zcr=n.hstack(zcr_a)

xcv=n.hstack(xcv_a)
ycv=n.hstack(ycv_a)
zcv=n.hstack(zcv_a)

xct=n.hstack(xct_a)
yct=n.hstack(yct_a)
zct=n.hstack(zct_a)

# define grid.
xi = n.logspace(-0.51, 0.51, 25) #n.arange(0.25, 3.2, 0.02)
yi = n.logspace(-0.51, 0.51, 25) #n.arange(0.25, 3.2, 0.02)
# grid the data.
cci = griddata((xcr, ycr), zcr, (xi[None,:], yi[:,None]), method='linear')
cvi = griddata((xcv, ycv), zcv, (xi[None,:], yi[:,None]), method='linear')

var_cr_array = n.array([ var_cr_L40(xcr[ii], ycr[ii]) for ii in range(len(xcr))])
cci_m = griddata((xcr, ycr), var_cr_array, (xi[None,:], yi[:,None]), method='linear')

var_total_array = n.array([ var_total_L40(xcv[ii], ycv[ii]) for ii in range(len(xcv))])
cvi_m = griddata((xcv, ycv), var_total_array, (xi[None,:], yi[:,None]), method='linear')

p.scatter(-n.log10(xcv), -n.log10(ycv), c=n.log10(var_total_array), s=10, edgecolors='none',vmin=-4,vmax=0)
#CS = p.contour(-n.log10(xi), -n.log10(yi), n.log10(cvi_m),15,linewidths=0.5,colors='k', vmin=-4, vmax=0)
#CS = p.contourf(-n.log10(xi), -n.log10(yi), n.log10(cvi_m),15,cmap=p.cm.jet, vmin=-4, vmax=0)
cb = p.colorbar(shrink=0.8)
cb.set_label("R($\sigma_1, \sigma_2$)")
p.xlabel(r'$log_{10}(\sigma_1^{-1})$')
p.ylabel(r'$log_{10}(\sigma_2^{-1})$')
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_model.png"))
p.clf()

p.scatter(-n.log10(xct), -n.log10(yct), c=n.log10(zct), s=10, edgecolors='none',vmin=-4,vmax=0)
#p.scatter(-n.log10(xcv), -n.log10(ycv), c=n.log10(zcv), s=10, edgecolors='none',vmin=-4,vmax=0)
#CS = p.contour(-n.log10(xi), -n.log10(yi), n.log10(cvi),15,linewidths=0.5,colors='k', vmin=-4, vmax=0)
#CS = p.contourf(-n.log10(xi), -n.log10(yi), n.log10(cvi),15,cmap=p.cm.jet, vmin=-4, vmax=0)
cb = p.colorbar(shrink=0.8)
cb.set_label("R($\sigma_1, \sigma_2$)")
p.xlabel(r'$log_{10}(\sigma_1^{-1})$')
p.ylabel(r'$log_{10}(\sigma_2^{-1})$')
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_data.png"))
p.clf()


p.scatter(-n.log10(xcr), -n.log10(ycr), c=n.log10(var_cr_array), s=10, edgecolors='none',vmin=-5,vmax=0)
#CS = p.contour(-n.log10(xi), -n.log10(yi), cci_m,15,linewidths=0.5,colors='k', vmin=0, vmax=1)
#CS = p.contourf(-n.log10(xi), -n.log10(yi), cci_m,15,cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label("R($\sigma_1, \sigma_2$)")
p.xlabel(r'$log_{10}(\sigma_1^{-1})$')
p.ylabel(r'$log_{10}(\sigma_2^{-1})$')
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"ccr_matrix_model.png"))
p.clf()

CS = p.contour(-n.log10(xi), -n.log10(yi), cci, 15, linewidths=0.5,colors='k', vmin=0, vmax=1)
CS = p.contourf(-n.log10(xi), -n.log10(yi), cci, 15, cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label("R($\sigma_1, \sigma_2$)")
p.xlabel(r'$log_{10}(\sigma_1^{-1})$')
p.ylabel(r'$log_{10}(\sigma_2^{-1})$')
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"ccr_matrix_data.png"))
p.clf()

sys.exit()


# contour the gridded data, plotting dots at the randomly spaced data points.
CS = p.contour(xi,yi,cci,15,linewidths=0.5,colors='k', vmin=0, vmax=1)
CS = p.contourf(xi,yi,cci,15,cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label("R($\sigma_1, \sigma_2$)")
p.xlabel(r'$\sigma_1$')
p.ylabel(r'$\sigma_2$')
p.plot(x,y,'k,')
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"griddata.png"))
p.clf()



# combines the covariance matrix
crC=n.zeros((79,79))
crItp=[]
cvC=n.zeros((79,79))

for ii in range(len(iis)):
	crC[ids[ii]:,ids[ii]:] = crs[ii][ids[ii]:,ids[ii]:]
	cvC[ids[ii]:,ids[ii]:] = cvs[ii][ids[ii]:,ids[ii]:]
	xx, yy = n.meshgrid(sigs[ii][ids[ii]:], sigs[ii][ids[ii]:])
	crItp.append(interp2d(xx, yy, crs[ii][ids[ii]:,ids[ii]:], kind='cubic', copy=True, bounds_error=False, fill_value=0.))

logmass = lgm[0]
fig = p.figure(0,(6,6))
mat = p.matshow(crC)
p.xticks(n.arange(0,len(logmass),5), logmass[n.arange(0,len(logmass),5)],rotation=45)
p.yticks(n.arange(0,len(logmass),5), logmass[n.arange(0,len(logmass),5)])
for iid in ids:
	p.axvline(iid-1, lw=1, color='k', ls='dashed')
	p.axhline(iid-1, lw=1, color='k', ls='dashed')

cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef")
p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.xlim((n.min(ids)-1, 75))
p.ylim((75, n.min(ids)-1))
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"mvir-cr-all.png"))
p.clf()

#### antidiagonal projection

def getXY(ipos = 40):
	delt = len(logmass) - ipos
	Ntot = len(logmass)
	antiD=[]
	xp=[]
	for ii in n.arange(0, Ntot-ipos, 1):
		antiD.append(crC[ipos+ii][ipos-ii])
		xp.append(logmass[ipos+ii])

	antiD = n.array(antiD)
	xp = n.array(xp)
	p.plot(xp, antiD, label='M='+str(logmass[ipos]))
	return xp, antiD

getXY(30)
getXY(35)
getXY(40)
getXY(45)
getXY(50)
getXY(55)
getXY(60)
getXY(65)
getXY(70)
p.ylabel("corrCoef")
p.xlabel('mass')
p.xlim((logmass[25], 16.))
p.ylim((-0.05, 1.05))
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"mvir-antiD.png"))
p.clf()


def getXYh(ipos = 40):
	delt = len(logmass) - ipos
	Ntot = len(logmass)
	antiD=[]
	xp=[]
	for ii in n.arange(0, Ntot-ipos, 1):
		antiD.append(crC[ipos+ii][ipos])
		xp.append(logmass[ipos+ii])

	antiD = n.array(antiD)
	xp = n.array(xp)
	p.plot(xp, antiD, label='M='+str(logmass[ipos]))
	return xp, antiD

getXYh(30)
getXYh(35)
getXYh(40)
getXYh(45)
getXYh(50)
getXYh(55)
getXYh(60)
getXYh(65)
getXYh(70)
p.ylabel("corrCoef")
p.xlabel('mass')
p.xlim((logmass[25], 16.))
p.ylim((-0.05, 1.05))
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"mvir-corrHoriZ.png"))
p.clf()

sss=n.logspace(-1,1,1000)
yyy = bias(sss)**(-0.5)
bmax = n.max(yyy)
bmin = n.min(yyy) 
amin = n.min(zi[n.isnan(zi)==False])
bfun = lambda xx : (bias(xx)**(-0.5)- bmin + amin)/(bmax-bmin+ amin)
xis, yis = n.meshgrid(xi,yi)


zm = bfun(xis) * bfun(yis)
CS = p.contour(xi,yi,zm,15,linewidths=0.5,colors='k', vmin=0, vmax=1)
CS = p.contourf(xi,yi,zm,15,cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label(r"$1/(\sqrt{b(\sigma_1)b(\sigma_2)})$")
p.xlabel(r'$\sigma$')
p.ylabel(r'$\sigma$')
p.plot(x,y,'k,')
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"griddata_model.png"))
p.clf()


CS = p.contour(xi,yi,zi-zm,15,linewidths=0.5,colors='k', vmin=0, vmax=1)
CS = p.contourf(xi,yi,zi-zm,15,cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label("residual CC - bias model")
p.xlabel(r'$\sigma$')
p.ylabel(r'$\sigma$')
p.plot(x,y,'k,')
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"griddata_resid.png"))
p.clf()

tr_val = griddata((x, y), z, (2.6, 2.6), method='linear')-bfun(2.6) * bfun(2.6)
bl_val = griddata((x, y), z, (0.5, 0.5), method='linear')-bfun(0.5) * bfun(0.5)

bnorm = lambda xx : 1-(((bias(xx) - bias(n.min(yi)) )/(bias(n.max(yi))-bias(n.min(yi))) )*(bl_val-tr_val)+tr_val )
 
def res_funs(aa, bb, yVal ): 
	return norm.pdf((bb-aa)/2, loc=0., scale=yVal/20. ) / norm.pdf(0. , loc=0., scale=yVal/20. ) * bnorm(yVal) 

zm2=n.zeros_like(yis)
for jj in range(len(yis)):
	zm2[jj] = res_funs(yis[jj], xis[jj], yi[jj]) 

CS = p.contour(xi,yi,zm2,15,linewidths=0.5,colors='k', vmin=0, vmax=1)
CS = p.contourf(xi,yi,zm2,15,cmap=p.cm.jet, vmin=0, vmax=1)
cb = p.colorbar(shrink=0.8)
cb.set_label("residual")
p.xlabel(r'$\sigma$')
p.ylabel(r'$\sigma$')
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"griddata_resid_model.png"))
p.clf()


for xid, xx in enumerate(xi):
	ztest = griddata((x, y), z, (xx, yi[yi<xx]), method='linear')
	p.plot(yi[yi<xx],ztest)#,label=str(n.round(xx,2)))

#p.plot(xi, bias(xi)**(-0.5)-(bmax-1),'k--', lw=2, label)
p.plot(xi, bias(xi)**(-0.5)-(bmax-1),'k--', lw=2, label=r"$1+b^{-0.5}_h-max(b^{-0.5}_h)$")
-derivative(bias,xi)
p.ylabel("corrCoef")
p.xlabel('sigma')
p.yscale('log')
p.ylim((0.001, 1.05))
gl = p.legend(loc=0,fontsize=14)
gl.set_frame_on(False)
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"sigma-corrHoriZ-log.png"))
p.clf()


for xid, xx in enumerate(xi):
	ztest = griddata((x, y), z, (xx, yi[yi<xx]), method='linear')
	p.plot(yi[yi<xx],ztest)#,label=str(n.round(xx,2)))

p.plot(xi, bias(xi)**(-0.5)-(bmax-1),'k--', lw=2, label=r"$1+b^{-0.5}_h-max(b^{-0.5}_h)$")
p.ylabel("corrCoef")
p.xlabel('sigma')
p.ylim((-0.05, 1.05))
gl = p.legend(loc=0,fontsize=14)
gl.set_frame_on(False)
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"sigma-corrHoriZ.png"))
p.clf()

"""
sss=n.logspace(-1,1,1000)
yyy = bias(sss)**(-0.5)
bmax = n.max(yyy)
p.plot(sss, )
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"bias-relation.png"))
p.clf()


sigmas = n.arange( 3.5, n.min(sigs), -0.02)
grid_x, grid_y = n.meshgrid(sigmas,sigmas)
#values = 
xx, yy = n.meshgrid(sigs[ii][ids[ii]:], sigs[ii][ids[ii]:])
#crItp.append(interp2d(xx, yy, crs[ii][ids[ii]:,ids[ii]:], 
out = crItp[0]( X,Y)#n.hstack((X)), n.hstack((Y)) )
#mat = out.reshape((len(sigmas),len(sigmas)))

crSS=n.zeros(( len(sigmas), len(sigmas) ))
sigMax = n.array([n.max(sigs[ii][ids[ii]:]) for ii in range(len(iis))])
for ii in range(len(iis)):
	print n.max(sigs[ii][ids[ii]:])
	sel=sigmas[(sigmas<n.max(sigs[ii][ids[ii]:]))]
	idx = n.argmax(sigmas<n.max(sigs[ii][ids[ii]:]))
	crSS[idx:, idx:] = crItp[ii](sel,sel)


fig = p.figure(0,(6,6))
mat = p.matshow(crSS)
p.xticks(n.arange(0,len(sigmas),5), n.round(sigmas[n.arange(0,len(sigmas),5)],3),rotation=45)
p.yticks(n.arange(0,len(sigmas),5), n.round(sigmas[n.arange(0,len(sigmas),5)],3))
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef")
p.xlabel(r'$\sigma$')
p.ylabel(r'$\sigma$')
p.xlim((130, 275))
p.ylim((130, 275))
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"sigma-cr-all.png"))
p.clf()

sys.exit()

#rebinned x 2
dataR = n.array([dt[2::2]+dt[1::2] for dt in data])
binsR = bins[1::2]
logmassR = ( binsR[1:]  + binsR[:-1] )/2.
NcountsR = dataR.sum(axis=0) 
okR= ( logmassR> logmp) & (NcountsR>2)

cvR = n.cov(dataR.T[okR])
crR = n.corrcoef(dataR.T[okR])
mmR = logmassR[okR]

mass2XR = interp1d(mmR, n.arange(len(mmR)))

fig = p.figure(0,(6,6))
mat = p.matshow(crR)
p.xticks(n.arange(0,len(mmR),5), mm[n.arange(0,len(mmR),5)],rotation=45)
p.yticks(n.arange(0,len(mmR),5), mm[n.arange(0,len(mmR),5)])
p.axvline(mass2XR(logmp+3), lw=2, color='k')
p.axhline(mass2XR(logmp+3), lw=2, color='k')
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef Mvir Hist Counts")
p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],"mvir-cr-2.png"))
p.clf()

#rebinning case
dataR = n.array([dt[2::2]+dt[1::2] for dt in data])
binsR = bins[1::2]
logmassR = ( binsR[1:]  + binsR[:-1] )/2.
NcountsR = dataR.sum(axis=0) 
okR= ( logmassR> logmp-0.5) & (NcountsR>2)

cvR = n.cov(dataR.T[okR])
crR = n.corrcoef(dataR.T[okR])
mmR = logmassR[okR]

mass2XR = interp1d(mmR, n.arange(len(mmR)))

fig = p.figure(0,(6,6))
mat = p.matshow(crR)
p.xticks(n.arange(0,len(mmR),5), mmR[n.arange(0,len(mmR),5)],rotation=45)
p.yticks(n.arange(0,len(mmR),5), mmR[n.arange(0,len(mmR),5)])
p.axvline(mass2XR(logmp+3), lw=2, color='k')
p.axhline(mass2XR(logmp+3), lw=2, color='k')
p.axvline(mass2XR(logmp+1), lw=2, color='k')
p.axhline(mass2XR(logmp+1), lw=2, color='k')
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef Mvir Hist Counts")
p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"mvir-cr-2_"+figName+".png"))
p.clf()



Ncounts = data.sum(axis=0) 
Nall = Ncounts / volume
ok= ( sigma> logmp-0.5) & (Ncounts>2)

index=n.arange(int(data.shape[0]))
n.random.shuffle( index )
Ntotal = int(data.shape[0])

dataS = n.array([n.sum(data[id:id+Ntotal/resamp:1], axis=0) for id in n.arange(0,Ntotal,Ntotal/resamp)])

cv = n.cov(data.T[ok])
cr = n.corrcoef(data.T[ok])
mm = logmass[ok]
sigma = sig[ok]
nu = nus[ok]

cvS = n.cov(dataS.T[ok])
crS = n.corrcoef(dataS.T[ok])

mass2X = interp1d(mm, n.arange(len(mm)))

fig = p.figure(0,(6,6))
mat = p.matshow(cr)
p.xticks(n.arange(0,len(mm),5), mm[n.arange(0,len(mm),5)],rotation=45)
p.yticks(n.arange(0,len(mm),5), mm[n.arange(0,len(mm),5)])
p.axvline(mass2X(logmp+3), lw=2, color='k')
p.axhline(mass2X(logmp+3), lw=2, color='k')
p.axvline(mass2X(logmp+1), lw=2, color='k')
p.axhline(mass2X(logmp+1), lw=2, color='k')
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef Mvir Counts "+figName)
p.xlabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[M_{vir}/(h^{-1}M_\odot)]$')
p.grid()
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"mvir-cr-0_"+figName+".png"))
p.clf()

fig = p.figure(0,(6,6))
mat = p.matshow(cv)
p.xticks(n.arange(0,len(nu),5), n.round(nu[n.arange(0,len(nu),5)],3),rotation=45)
p.yticks(n.arange(0,len(nu),5), n.round(nu[n.arange(0,len(nu),5)],3))
#p.axvline(mass2X(logmp+3), lw=2, color='k')
#p.axhline(mass2X(logmp+3), lw=2, color='k')
#p.axvline(mass2X(logmp+1), lw=2, color='k')
#p.axhline(mass2X(logmp+1), lw=2, color='k')
cb = p.colorbar(shrink=0.8)
cb.set_label("corrCoef Mvir Counts "+figName)
p.xlabel(r'$\nu$')
p.ylabel(r'$\nu$')
p.grid()

"""
