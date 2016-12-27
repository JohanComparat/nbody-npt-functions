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


fun = lambda lg_X, lg_A, lg_X0, lg_alpha, lg_beta : n.log10( 10**lg_A * (10**lg_X/10**lg_X0)**(-10**lg_beta) * n.e**(- (10**lg_X/10**lg_X0)**(10**lg_alpha) ) )
funG = lambda lg_X, lg_z, ps  : fun( lg_X, ps[0], ps[1], ps[2], ps[3] ) #
	
NDecimal = 3
dir='..'
qty = 'M200c'
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")
dir_25N = join(dir,"MD_2.5GpcNW")
dir_40N = join(dir,"MD_4GpcNW")

data = fits.open( join(dir, qty, "MD_M200c_summary.fits") )[1].data


errorLog = 0.03
NminCount = 0
Npmin = 3
limits_04 = [Npmin*9.63 * 10**7, 5e12]
limits_10 = [Npmin*1.51 * 10**9., 5e13]
limits_25 = [Npmin*2.359 * 10**10., 5e14]
limits_40 = [Npmin* 9.6 * 10**10. , 5e15]
MPART = n.array([9.63 * 10**7, 1.51 * 10**9, 2.359 * 10**10, 9.6 * 10**10])
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin = -0.001
zmax = 0.001

def GetLim(xx , ratio, xTr , lims):
	lowX = (xx<xTr)
	highX = (xx>xTr)
	#X_up = n.array([ n.min(xx[(highX)&(ratio<lim)]) for lim in lims ])
	X_low = n.array([ n.max(xx[(lowX)&(ratio<lim)]) for lim in lims ])
	return X_low # X_up,


def getCompleteness(qty = 'M200c', cos = "cen", zmin = -0.01, zmax = 0.1):
	"""
	Plots the data to be used in the fits later in the analysis.
	"""
	# gets the best fitting parameters at redshift 0 :
	f=open(join(dir,qty,"M200c-"+cos+"-cumulative-function-z0-params.pkl"), 'r')
	res = cPickle.load(f)
	f.close()
	pOpt = res.x
	cov = res.direc
	
	# redshift selection
	zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
	# minimum number counts selection
	nSel = (data['dN_counts_'+cos]>NminCount)
	
	#start the figure
	p.figure(1,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	
	# mass selection and plots for each box :
	if  cos == "cen":
		mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0]))) 
		ok = (zSel) & (mSel) & (nSel)
		lg_M200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
		lg_MF_c = n.log10(data["dNdVdlnM_"+cos+"_c"][ok])
		lg_1pz = n.log10(1+ data["redshift"][ok])
		
		X_low_04 = GetLim(xx = lg_M200c, ratio = 10**(lg_MF_c - funG(lg_M200c,lg_1pz, pOpt)), xTr = 2*limits_04[0], lims = [0.8, 0.9, 0.95, 0.97])
		p.plot(lg_M200c[::3], 10**(lg_MF_c[::3] - funG(lg_M200c[::3],lg_1pz[::3], pOpt)), marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)

		mSel = ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0]))) 
		ok = (zSel) & (mSel) & (nSel)
		lg_M200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
		lg_MF_c = n.log10(data["dNdVdlnM_"+cos+"_c"][ok])
		lg_1pz = n.log10(1+ data["redshift"][ok])
		
		X_low_10 = GetLim(xx = lg_M200c, ratio = 10**(lg_MF_c - funG(lg_M200c,lg_1pz, pOpt)), xTr = 2*limits_10[0], lims = [0.8, 0.9, 0.95, 0.97])
		p.plot(lg_M200c[::3], 10**(lg_MF_c[::3] - funG(lg_M200c[::3],lg_1pz[::3], pOpt)), marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)

		
		mSel = ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>n.log10(limits_25[0]))) 
		ok = (zSel) & (mSel) & (nSel)
		lg_M200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
		lg_MF_c = n.log10(data["dNdVdlnM_"+cos+"_c"][ok])
		lg_1pz = n.log10(1+ data["redshift"][ok])
		
		X_low_25 = GetLim(xx = lg_M200c, ratio = 10**(lg_MF_c - funG(lg_M200c,lg_1pz, pOpt)), xTr = 2*limits_25[0], lims = [0.8, 0.9, 0.95, 0.97])
		p.plot(lg_M200c[::3], 10**(lg_MF_c[::3] - funG(lg_M200c[::3],lg_1pz[::3], pOpt)), marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)


		mSel = ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>n.log10(limits_40[0]))) 
		ok = (zSel) & (mSel) & (nSel)
		lg_M200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
		lg_MF_c = n.log10(data["dNdVdlnM_"+cos+"_c"][ok])
		lg_1pz = n.log10(1+ data["redshift"][ok])
		
		X_low_40 = GetLim(xx = lg_M200c, ratio = 10**(lg_MF_c - funG(lg_M200c,lg_1pz, pOpt)), xTr = 2*limits_40[0], lims = [0.8, 0.9, 0.95, 0.97])
		p.plot(lg_M200c, 10**(lg_MF_c - funG(lg_M200c,lg_1pz, pOpt)), marker ='+', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)

		p.axvline(n.log10(MPART[0]*100.), color='r')
		p.axvline(n.log10(MPART[1]*100.), color='c')
		p.axvline(n.log10(MPART[2]*100.), color='m')
		p.axvline(n.log10(MPART[3]*100.), color='b')
		p.xlim((9.5,16))
		
		lims = n.array([X_low_04, X_low_10, X_low_25, X_low_40]) # [X_up_04, X_up_10, X_up_25, X_up_40],

		
		
		
	if  cos == "sat":
		mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0]))) 
		ok = (zSel) & (mSel) & (nSel)
		lg_M200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
		lg_MF_c = n.log10(data["dNdVdlnM_"+cos+"_c"][ok])
		lg_1pz = n.log10(1+ data["redshift"][ok])
		
		X_low_04 = GetLim(xx = lg_M200c, ratio = 10**(lg_MF_c - funG(lg_M200c,lg_1pz, pOpt)), xTr = 2*limits_04[0], lims = [0.8, 0.9, 0.95, 0.97])
		p.plot(lg_M200c[::3], 10**(lg_MF_c[::3] - funG(lg_M200c[::3],lg_1pz[::3], pOpt)), marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)

		mSel = ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0]))) 
		ok = (zSel) & (mSel) & (nSel)
		lg_M200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
		lg_MF_c = n.log10(data["dNdVdlnM_"+cos+"_c"][ok])
		lg_1pz = n.log10(1+ data["redshift"][ok])
		
		X_low_10 = GetLim(xx = lg_M200c, ratio = 10**(lg_MF_c - funG(lg_M200c,lg_1pz, pOpt)), xTr = 2*limits_10[0], lims = [0.8, 0.9, 0.95, 0.97])
		p.plot(lg_M200c[::3], 10**(lg_MF_c[::3] - funG(lg_M200c[::3],lg_1pz[::3], pOpt)), marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)

		p.axvline(n.log10(MPART[0]*100.), color='r')
		p.axvline(n.log10(MPART[1]*100.), color='c')
		p.xlim((9.5,16))
		
		lims = n.array([X_low_04, X_low_10]) # [X_up_04, X_up_10, X_up_25, X_up_40],

	p.axhline(1)
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r' n('+cos+',>M) data / model') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=4,fontsize=12)
	gl.set_frame_on(False)

	p.ylim((.9,1.1))
	#p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0[0])+" "+str(vcut0[0])+" "+str(a0[0])+" "+str(b0[0]))
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-completeness.png"))
	p.clf()
	
	return lims
	
limC = getCompleteness(qty = 'M200c', cos = "cen", zmin = -0.01, zmax = 0.01)
limS = getCompleteness(qty = 'M200c', cos = "sat", zmin = -0.01, zmax = 0.01)

f=open( join( dir, qty, "completeness-" + qty + "-npart-z0.txt" ), 'w' )
f.write( " central \n")

for ii, el in enumerate(limC):
	f.write( names[ii]+" & " +str(n.round(n.log10(el[0]),NDecimal))+ " ("+str(int(el[0]/MPART[ii]))+ ") & " + str(n.round(n.log10(el[1]),NDecimal))+ " ("+str(int(el[1]/MPART[ii]))+ ") & " + str(n.round(n.log10(el[2]),NDecimal))+ " ("+str(int(el[2]/MPART[ii]))+") & " + str(n.round(n.log10(el[3]),NDecimal))+ " ("+str(int(el[3]/MPART[ii]))+ ") \\\\ \n")

f.write( " sat \n")

for ii, el in enumerate(limS):
	f.write( names[ii]+" & " +str(n.round(n.log10(el[0]),NDecimal)) + " ("+str(int(el[0]/MPART[ii]))+ ") & " + str(n.round(n.log10(el[1]),NDecimal)) + " ("+str(int(el[1]/MPART[ii]))+ ") & " + str(n.round(n.log10(el[2]),NDecimal)) + " ("+str(int(el[2]/MPART[ii]))+ ") & " + str(n.round(n.log10(el[3]),NDecimal)) + " ("+str(int(el[3]/MPART[ii]))+ ")\\\\ \n")

f.close()

