import glob
import sys
import cPickle
from os.path import join
import numpy as n
import astropy.io.fits as fits
import os

import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p
from scipy.optimize import minimize
from scipy.optimize import curve_fit

fun = lambda lg_X, lg_A, lg_X0, lg_alpha, lg_beta : n.log10( 10**lg_A * (10**lg_X/10**lg_X0)**(-10**lg_beta) * n.e**(- (10**lg_X/10**lg_X0)**(10**lg_alpha) ) )


dir='..'
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")
dir_25N = join(dir,"MD_2.5GpcNW")
dir_40N = join(dir,"MD_4GpcNW")

data = fits.open( join("..", "M200c", "MD_M200c_summary.fits") )[1].data


errorLog = 0.03
NminCount = 10
Npmin = 300
limits_04 = [Npmin*9.63 * 10**7, 5e12]
limits_10 = [Npmin*1.51 * 10**9., 5e13]
limits_25 = [Npmin*2.359 * 10**10., 5e14]
limits_40 = [Npmin* 9.6 * 10**10. , 5e15]
MPART = n.array([9.63 * 10**7, 1.51 * 10**9, 2.359 * 10**10, 9.6 * 10**10])
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin = -0.01
zmax = 0.01

def fitData(qty = 'M200c', cos = "cen", zmin = -0.01, zmax = 0.1, p0 = [-4., 13.0, -0.3, -0.04]):
	"""
	Plots the data to be used in the fits later in the analysis.
	"""
	# redshift selection
	zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
	# mass selection
	if  cos == "cen":
		mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0]))) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0]))) | ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>n.log10(limits_25[0]))) | ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>n.log10(limits_40[0]))) 
	if  cos == "sat":
		mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0]))) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])))  #| ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>n.log10(limits_25[0]))) | ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>n.log10(limits_40[0]))) 
	# minimum number counts selection
	nSel = (data['dN_counts_'+cos]>NminCount)
	# altogether
	ok = (zSel) & (mSel) & (nSel)
	# now the plot
	lg_M200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
	#print len(lg_M200c), lg_M200c
	lg_MF_c = n.log10(data["dNdVdlnM_"+cos+"_c"][ok])
	#print lg_MF_c
	lg_1pz = n.log10(1+ data["redshift"][ok])
	#print lg_1pz
	
	funG = lambda lg_X, lg_z, ps  : fun( lg_X, ps[0], ps[1], ps[2], ps[3] ) #
	chi2fun = lambda ps : n.sum( (funG(lg_M200c, lg_1pz, ps) - lg_MF_c)**2. / (errorLog)**2. )/(len(lg_MF_c) - len(ps))
	res = minimize(chi2fun, p0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
	pOpt = res.x
	cov = res.direc
	chi2perpoint = lambda ps : (funG(lg_M200c, lg_1pz, ps) - lg_MF_c)**2. / (errorLog)**2. 
	chi2pp = chi2perpoint(pOpt)
	print pOpt, cov

	lg_M200c_model = n.arange(n.min(lg_M200c),n.max(lg_M200c),0.1)
	X,Y = n.meshgrid(lg_M200c_model, n.arange(zmin, zmax+0.025,0.025))
	Z = funG(X,n.log10(1+Y),pOpt)
	n.savetxt(join(dir,qty,"M200c-"+cos+"-cumulative-function-z0-model-pts.txt"),n.transpose([n.hstack((X)), n.hstack((Y)), n.hstack((Z))]) )
	
	f=open(join(dir,qty,"M200c-"+cos+"-cumulative-function-z0-params.pkl"), 'w')
	cPickle.dump(res, f)
	f.close()
	
	X,Y,Z = n.loadtxt(join(dir,qty,"M200c-"+cos+"-cumulative-function-z0-model-pts.txt"), unpack=True)
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(X, Z, c=Y, s=5, marker='o',label="model", rasterized=True)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r'log  n(>M)') 
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((-8, 1))
	p.xlim((9.5,16))
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-model.png"))
	p.clf()

	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(lg_M200c, lg_MF_c, c=chi2pp, s=5, marker='o',label="chi2", rasterized=True)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("chi2 per point")
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r'log  n(>M)') 
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((-8, 1))
	p.xlim((9.5,16))
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-chi2PP.png"))
	p.clf()
	
fitData(qty = 'M200c', cos = "cen", zmin = -0.01, zmax = 0.1, p0 = [-4, 13.5, -0.2, -0.1])
fitData(qty = 'M200c', cos = "sat", zmin = -0.01, zmax = 0.1, p0 = [-4., 12.8, -0.3, -0.03])
