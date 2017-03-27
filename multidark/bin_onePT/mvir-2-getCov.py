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


qty = 'mvir'
dir = join(os.environ['MVIR_DIR'])
plotDir = "covariance"
# loads summary file
dataMF = fits.open( join(dir, qty+"_summary.fits"))[1].data
zzero = (dataMF['redshift']==0) & (dataMF['log_mvir']>3+dataMF['logMpart']) & (dataMF['dN_counts_cen'] > 10 )
dlnSigM = abs(n.log(dataMF['log_mvir_max']-dataMF['log_mvir_min'])*dataMF['dlnsigmaMdlnM'])

#lib.covariance_factor
#lib.f_BH(sigma, 0.333, 0.788, 0.807, 1.795)
biasHMF = lambda sigma : lib.b_BH(sigma, a=0.90343, p=0.64031, q=1.69561)
biasB = lambda sigma : lib.b_BH(sigma, 0.76966, 0.49524, 1.34957)

bias = biasB
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
version = 'v4'
# measurement files
fileC = n.array(glob.glob( join(os.environ['MD_DIR'], "MD_*Gpc*",  version, qty,"out_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MD_DIR'], "MD_*Gpc*", version, qty,"out_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MD_DIR'], "MD_*Gpc*", version, qty,"out_*_Satellite_JKresampling.pkl")))
fileC.sort()
fileB.sort()
fileS.sort()

# redshift 0 data
iis = [8, 13, 27, 21, -30, 30]#[-1, -2, -4, -9, -22, 3]
#n.arange(6) # [8, 13, 29, 31, 60, -10]
print fileC

lbox = n.array([40., 100., 250., 400. ])
vols = lbox**3.
snFactor = n.array([4.7, 5.4, 8.0, 10.0*0.95 ])
svFactor = 1.5*n.ones(4) #n.array([1.56, 1.61, 1.65, 1.414 ])
binW = 0.025 
	
p.figure(1,(6,6))
p.axes([0.17, 0.17, 0.75, 0.78])
xx = n.arange(1.,3., 0.1)
coefs=n.polyfit(n.log10(lbox), snFactor,1,w=10*n.ones_like(lbox))
p.plot(n.log10(lbox), snFactor, 'bo', label='MD covariance')
p.plot(xx,n.polyval(coefs,xx),'b--', label=r"$"+str(n.round(coefs[1],2))+"+"+str(n.round(coefs[0],2))+r'\, \log_{10}(L_{box})$')
#coefs=n.polyfit(n.log10(lbox), svFactor,1,w=10*n.ones_like(lbox))
#p.plot(n.log10(lbox), svFactor, 'ro', label='sv')
#p.plot(xx,n.polyval(coefs,xx),'r--', label=str(n.round(coefs,2)))
p.xlabel(r'$\log_{10}(L_{box})$')
p.ylabel(r'$Q$')#,\; f_{sv}
p.grid()
gl=p.legend(loc=0, frameon=False)
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_sn_sv_lbox.png"))
p.clf()

p.figure(1,(6,6))
p.axes([0.17, 0.17, 0.75, 0.78])
p.plot(n.log10(vols), snFactor, 'bo', label='sn')
#p.plot(n.log10(vols), svFactor, 'ro', label='sv')
p.xlabel(r'$\log_{10}(V_{box})$')
p.ylabel(r'$Q$') # ,\; f_{sv}
p.grid()
gl=p.legend(loc=0, frameon=False)
p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_sn_sv_vol.png"))
p.clf()

#sys.exit()

def plot_COV(fileC, binFile):
	"""
	Reads the data and extract mass function covariance matrix
	
	First gets the basic info about the data
	
	Then creates a dimension-less mass function matrix (1000 mass functions)
	
	"""
	# hig mass end M25 and M40
	#low mass end M04 and M10
	#snFactor = n.array([4.8/1.1,    4.2/0.95,  6.6, 9.1, 12])
	#svFactor = n.array([2./1.1, 3.5/0.95 , 3.1,  5,    4.5])
	#binW = 0.125 
	#plotDir = "covariance_lt14"

	ddd, boxName = lib.convert_pkl_massFunction_covarianceMatrix(fileC, binFile, qty='mvir', delta_wrt='mean')
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
	p.axes([0.17, 0.17, 0.75, 0.78])
	p.title(nickname[boxName])#+" total")
	tp=(xcv>ycv)# (xcv!=ycv)
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=n.log10(ctotal[tp]), s=60, edgecolors='none',vmin=-2,vmax=2., marker='s', label='data')
	tp=(xcv<ycv)# (xcv!=ycv)
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=n.log10(model[boxName](xcv, ycv)[tp]), s=60, edgecolors='none', vmin=-2,vmax=2, marker='o', label='model')
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
	tp=(xcv>=ycv)
	residuals = n.log10(ctotal[tp]/model[boxName](xcv, ycv)[tp])
	p.scatter(-n.log10(xcv[tp]), -n.log10(ycv[tp]), c=residuals, s=60, edgecolors='none', vmin=-0.5,vmax=0.5, marker='s')
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
	
	
	zSel = (zzero) & (dataMF['boxName'] ==boxName)
		
	p.figure(2,(6,6))
	p.axes([0.17, 0.17, 0.78, 0.78])
	p.title(nickname[boxName])
	y=model_sn[boxName](xcv, ycv)
	p.plot(-n.log10(xcv.diagonal()), n.log10(y.diagonal()), 'k-.', label='shot noise')
	y=model_sv[boxName](xcv, ycv)
	p.plot(-n.log10(xcv.diagonal()), n.log10(y.diagonal()), 'k--', label='sample variance')
	y=model[boxName](xcv, ycv)
	p.plot(-n.log10(xcv.diagonal()), n.log10(y.diagonal()), 'k', label='$C_{model}$')

	p.plot(-n.log10(xcv.diagonal()), n.log10(ctotal.diagonal()), 'b+', label=r'$C$')
	p.plot(-n.log10(dataMF['sigmaM'][zSel]), n.log10(dataMF["std90_pc_cen"][zSel]), 'rx', label=r'$C^{JK}$')

	#p.xlim((-0.7, 0.6))
	p.grid()
	#p.ylim((-0.7, 0.6))
	p.xlabel(r'$log_{10}(\sigma^{-1})$')
	p.ylabel(r'$log_{10}(C,\; C_{JK})$')
	gl=p.legend(loc=0, frameon=False)
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_"+boxName+"_diagonal.png"))
	p.clf()
	
	p.figure(2,(6,6))
	p.axes([0.17, 0.17, 0.78, 0.78])
	p.title(nickname[boxName])
	y=model[boxName](xcv, ycv).diagonal() / ctotal.diagonal()
	p.plot(-n.log10(xcv.diagonal()), y, 'k', label='C')
	interpData = interp1d(-n.log10(dataMF['sigmaM'][zSel]), dataMF["std90_pc_cen"][zSel])
	y2= model[boxName](xcv, ycv).diagonal()/ interpData(-n.log10(xcv.diagonal()))
	p.plot(-n.log10(xcv.diagonal()), y2, 'b', label='CJK')
	
	if boxName=="MD_0.4Gpc" or boxName=="MD_1Gpc" :
		sel = (-n.log10(xcv.diagonal())<-0.2)
	if boxName=="MD_2.5Gpc" or boxName=="MD_2.5GpcNW" :
		sel = (-n.log10(xcv.diagonal())<0.1)
	if boxName=="MD_4Gpc" or boxName=="MD_4GpcNW" :
		sel = (-n.log10(xcv.diagonal())<0.13)
		
	y3= model[boxName](xcv, ycv).diagonal()[sel]/ interpData(-n.log10(xcv.diagonal()))[sel]
	p.axhline(n.median(y3), label=str(n.round(n.median(y3),2)))
	#p.xlim((-0.7, 0.6))
	p.grid()
	p.ylim((0.1, 2.5))
	p.xlabel(r'$log_{10}(\sigma^{-1})$')
	p.ylabel(r'diagonal $C_{model}/C_{data}$')
	gl=p.legend(loc=0, frameon=False)
	p.savefig(join(os.environ['MVIR_DIR'],plotDir,"covariance_matrix_"+boxName+"_diagonal_residual.png"))
	p.clf()
	
	
	return f_mean, f_matrix, count_matrix, sigma, mass, residuals

def compCov():
	dout = []
	plotDir = "covariance"
	MDnames= n.array(['M04', 'M10', 'M25','M25n','M40','M40n'])# , 'D80'])
	for ii in iis:
		dout.append( plot_COV(fileC[ii], fileB[ii] ))

	"""
	fileC = n.array(glob.glob( join(os.environ['DS_DIR'], version, qty,"ds*_Central_JKresampling.pkl")))
	fileB = n.array(glob.glob( join( os.environ['DS_DIR'], version, qty,"ds*_"+qty+"_JKresampling.bins")))
	fileS = n.array(glob.glob( join( os.environ['DS_DIR'], version, qty,"ds*_Satellite_JKresampling.pkl")))
	plot_COV(fileC[0], fileB[0])
	ii=0
	dout.append( plot_COV(fileC[ii], fileB[ii]) )
	"""
	f_mean, f_matrix, count_matrix, sigma, mass, residuals = n.transpose(dout)
	

	#MDnames= n.array(['M04', 'M10', 'M25','M25n','M40','M40n', 'D80'])

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
	for nm, xm in  zip(MDnames, xmax):
		print nm, 10**xm

compCov()

