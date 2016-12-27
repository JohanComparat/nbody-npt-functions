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

dir='..'
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")
dir_25N = join(dir,"MD_2.5GpcNW")
dir_40N = join(dir,"MD_4GpcNW")

data = fits.open( join("..", "M200c", "MD_M200c_summary.fits") )[1].data

NDecimal = 3
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

def fitDataAll(qty = 'M200c', cos = "cen", zmin = -0.01, zmax = 2.3, p0 =  n.array([ 0.65, 12.2, -23.3, -0.9, -9.8, 14.9, -0.2, 1.23, -6.7, -11.6, 0.03, -0.33, 1.3 ])):
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
	# data selection
	ok = (zSel) & (mSel) & (nSel)
	# axis definition
	lg_M200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
	lg_MF_c = n.log10(data["dNdVdlnM_"+cos+"_c"][ok])
	lg_1pz = n.log10(1+ data["redshift"][ok])
	# fitting function definition
	# loads redshift 0 best parameters
	f=open(join(dir,qty,"M200c-"+cos+"-cumulative-function-z0-params.pkl"), 'r')
	res = cPickle.load(f)
	f.close()
	pInit = res.x
	print pInit
	# create redshift varying functions for the parameters
	A_pr = lambda lz, A1, A2, A3 : 10**(pInit[0] + A1 * lz + A2 * lz**2. + A3 * lz**3.)
	M_pr = lambda lz, m1, m2, m3 : 10**(pInit[1] + m1 * lz + m2 *lz**2. + m3 * lz**3.)
	a_pr = lambda lz, a1, a2, a3, a4 : 10**(pInit[2] + a1 * lz + a2 *lz**2. + a3 *lz**3. + a4 *lz**4.)
	b_pr = lambda lz, b1, b2, b3 : -10**(pInit[3] + b1 * lz + b2 *lz**2.+ b3 *lz**3.)
	# generalized fitting function
	vfG = lambda lg_v, lg_z, ps  : n.log10( A_pr(lg_z, ps[0], ps[1], ps[2]) * (10**lg_v/M_pr(lg_z, ps[3], ps[4], ps[5]))**b_pr(lg_z, ps[10], ps[11], ps[12])  * n.e**(- (10**lg_v/M_pr(lg_z, ps[3], ps[4], ps[5]))**a_pr (lg_z, ps[6], ps[7], ps[8], ps[9]) ) )

	# defines chi2
	chi2fun = lambda ps : n.sum( (vfG(lg_M200c, lg_1pz, ps) - lg_MF_c)**2. / (errorLog)**2. )/(len(lg_MF_c) - len(ps))
	# fits the parameters 
	res = minimize(chi2fun, p0, method='Powell',options={'xtol': 1e-8, 'disp': True, 'maxiter' : 5000000000000})
	pOpt = res.x
	cov = res.direc
	chi2perpoint = lambda ps : (vfG(lg_M200c, lg_1pz, ps) - lg_MF_c)**2. / (errorLog)**2. 
	chi2pp = chi2perpoint(pOpt)
	print pOpt

	lg_M200c_model = n.arange(n.min(lg_M200c),n.max(lg_M200c)+0.1,0.1)
	X,Y = n.meshgrid(lg_M200c_model, n.arange(zmin, zmax+0.02,0.02))
	Z = vfG(X,n.log10(1+Y),pOpt)
	outPointFile = join(dir,qty,"M200c-"+cos+"-cumulative-function-z0-model-pts.txt")
	n.savetxt(outPointFile ,n.transpose([n.hstack((X)), n.hstack((Y)), n.hstack((Z))]) )
	
	f=open(join(dir,qty,"M200c-"+cos+"-cumulative-function-zAll-params.pkl"), 'w')
	cPickle.dump(res, f)
	f.close()
	
	X,Y,Z = n.loadtxt(outPointFile, unpack=True)
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
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-model-zAll.png"))
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
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-chi2PP-zAll.png"))
	p.clf()
	# saves table of parameters fitted 
	f=open(join(dir,qty, "latex-parameters-"+qty+"-"+cos+".txt"),'w')
	f.write( "$A(z)$  & " + str(n.round(pInit[0],NDecimal))+ " & "+str(n.round(pOpt[0],NDecimal))+" & "+str(n.round(pOpt[1],NDecimal))+ " & "+str(n.round(pOpt[2],NDecimal))+ "\\\\ \n")
	f.write( "$M_{cut}(z)$ & " + str(n.round(pInit[1],NDecimal))+ " & "+str(n.round(pOpt[3],NDecimal))+" & "+str(n.round(pOpt[4],NDecimal))+ " & "+str(n.round(pOpt[5],NDecimal))+ "\\\\ \n")
	f.write( "$\alpha(z)$ & " + str(n.round(pInit[2],NDecimal))+ " & "+str(n.round(pOpt[6],NDecimal))+" & "+str(n.round(pOpt[7],NDecimal))+ " & "+str(n.round(pOpt[8],NDecimal))+ " & "+str(n.round(pOpt[9],NDecimal))+ "\\\\ \n")
	f.write( "$\beta(z)$ & " + str(n.round(pInit[3],NDecimal))+ " & "+str(n.round(pOpt[10],NDecimal))+" & "+str(n.round(pOpt[11],NDecimal))+ " & "+str(n.round(pOpt[12],NDecimal))+ "\\\\ \n")
	f.close()
	
	

print "centrals"	
fitDataAll(qty = 'M200c', cos = "cen", zmin = -0.01, zmax = 2.3, p0 =n.array([ 0.5, 13.1, -20.8, -0.89, -10.4, 13.6, -0.2, 0.84, -4.1, 6.5, 0.11, -0.82, 1.77 ]) )
print "satellites"
fitDataAll(qty = 'M200c', cos = "sat", zmin = -0.01, zmax = 2.3, p0 = n.array([ 0.8, 12., -22., -0.88, -9.8, 15., -0.28, 0.6, -0.3, 0.86, 0.15, -0.9, 1.8]))
