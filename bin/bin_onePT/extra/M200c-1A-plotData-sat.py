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

from scipy.interpolate import interp1d
from scipy.misc import derivative

sigma = n.arange(0.05,10,0.05)
f = lambda sigma, A, a, b, c : A*( (sigma/b)**(-a) + 1 )*n.e**(-c/sigma**2.)
ft08Y = f(sigma, 0.224, 1.67, 1.80, 1.48) 
ft08X = n.log10(1./sigma)

# mass and mass function from Klypin 16
k16x, k16y=n.loadtxt(join("..", "M200c","previousData","Klypin_2016_MF_z0.txt"),  unpack=True)

dir='..'
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")
dir_25N = join(dir,"MD_2.5GpcNW")
dir_40N = join(dir,"MD_4GpcNW")

data = fits.open( join("..", "M200c", "MD_M200c_summary.fits") )[1].data

boxes = set(data['boxName'])
mk = {"MD_0.4Gpc": '1', "MD_1Gpc": 'x' ,"MD_2.5Gpc": '2',"MD_4Gpc": '3' ,"MD_2.5GpcNW": '4', "MD_4GpcNW": '+'}

NminCount = 10
Npmin = 300
limits_04 = [Npmin*9.63 * 10**7, 5e12]
limits_10 = [Npmin*1.51 * 10**9., 5e13]
limits_25 = [Npmin*2.359 * 10**10., 5e14]
limits_40 = [Npmin* 9.6 * 10**10. , 5e15]
MPART = n.array([9.63 * 10**7, 1.51 * 10**9, 2.359 * 10**10, 9.6 * 10**10])
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin = -0.01
zmax = 0.001

def plotData(qty = 'M200c', cos = "cen", zmin = -0.01, zmax = 2.3):
	"""
	Plots the data to be used in the fits later in the analysis.
	"""
	# redshift selection
	zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
	# mass selection
	mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0])) &(data["log_"+qty+"_max"]<n.log10(limits_04[1]))) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])) &(data["log_"+qty+"_max"]<n.log10(limits_10[1]))) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>n.log10(limits_25[0])) &(data["log_"+qty+"_max"]<n.log10(limits_25[1]))) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>n.log10(limits_40[0]))&(data["log_"+qty+"_max"]<n.log10(limits_40[1]))) 
	#if  cos == "sat":
	#mSel =  ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0])) &(data["log_"+qty+"_max"]<n.log10(limits_04[1]))) |  ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])) &(data["log_"+qty+"_max"]<n.log10(limits_10[1])))
	# minimum number counts selection
	nSel = (data['dN_counts_'+cos]>NminCount)
	# altogether
	ok = (zSel) & (mSel) & (nSel)
	
	# x coordinates
	logsigM1 = n.log10(1./data['sigmaM'][ok])#
	print n.min(logsigM1), n.max(logsigM1)
	log_m200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
	m200c = 10**log_m200c
	# mean density array normalization
	rhom = cosmo.critical_density(data["redshift"][ok]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][ok])/(100*uu.km/(uu.Mpc*uu.s)))**1.
	# y coordinates
	log_MF = n.log10( m200c * data["dNdVdlnM_"+cos][ok]/ rhom.value )
	log_MF_c = n.log10(  data["dNdVdlnM_"+cos+"_c"][ok])
	log_f =  n.log10(m200c * data["dNdVdlnM_"+cos][ok]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"][ok]))
	log_f_c =  n.log10(m200c * data["dNdVdlnM_"+cos+"_c"][ok]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"][ok]))
	
	# now the plots
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	for box in boxes:
		bb = (data["boxName"][ok]==box)
		#print box, mk[box]
		sc1=p.scatter(logsigM1[bb], log_MF[bb], c=data["redshift"][ok][bb], s=20, marker=mk[box],label=box, rasterized=True, vmin=zmin, vmax = zmax, edgecolors='face')
	
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'$ln(\sigma^{-1})$')
	p.ylabel(r'log$_{10} (M^2/\rho_m) dn(M)/dM$') 
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.xlim((-0.7,0.6))
	p.ylim((-4.5,-2))
	#p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-differential-function-data-xSigma.png"))
	p.clf()

	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	for box in boxes:
		bb = (data["boxName"][ok]==box)
		#print box, mk[box]
		sc1=p.scatter(log_m200c[bb], log_MF[bb], c=data["redshift"][ok][bb], s=20, marker=mk[box],label=box, rasterized=True, vmin=zmin, vmax = zmax, edgecolors='face')
	
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r'log$_{10} [(M^2/\rho_m) dn(M)/dM]$') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	p.xlim((9.5,16))
	p.ylim((-4.5,-2))
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-differential-function-data.png"))
	p.clf()
	"""
	MD04=(data["boxLength"]==400.)
	MD10=(data["boxLength"]==1000.)
	MD25=(data["boxLength"]==2500.)
	MD40=(data["boxLength"]==4000.)
	log_m200c_04 = (data["log_"+qty+"_min"][ok & MD04]+data["log_"+qty+"_max"][ok & MD04])/2.
	error_04 = data['std90_pc_cen'][ok & MD04]
	
	log_m200c_10 = (data["log_"+qty+"_min"][ok & MD10]+data["log_"+qty+"_max"][ok & MD10])/2.
	error_10 = data['std90_pc_cen'][ok & MD10]
	
	log_m200c_25 = (data["log_"+qty+"_min"][ok & MD25]+data["log_"+qty+"_max"][ok & MD25])/2.
	error_25 = data['std90_pc_cen'][ok & MD25]
	
	log_m200c_40 = (data["log_"+qty+"_min"][ok & MD40]+data["log_"+qty+"_max"][ok & MD40])/2.
	error_40 = data['std90_pc_cen'][ok & MD40]
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_m200c_04, 100*error_04, c=data["redshift"][ok & MD04], s=5, marker='o',label="MD04", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	p.xlim((9.5,16))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-data04-uncertainty.png"))
	p.clf()
	
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_m200c_10, 100*error_10, c=data["redshift"][ok & MD10], s=5, marker='o',label="MD10", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	p.xlim((9.5,16))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-data10-uncertainty.png"))
	p.clf()
	
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_m200c_25, 100*error_25, c=data["redshift"][ok & MD25], s=5, marker='o',label="MD25", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	p.xlim((9.5,16))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-data25-uncertainty.png"))
	p.clf()
	
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_m200c_40, 100*error_40, c=data["redshift"][ok & MD40], s=5, marker='o',label="MD40", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r'Jackknife resampling relative error [%]') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	p.ylim((2e-2,30))
	p.xlim((9.5,16))
	p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-data40-uncertainty.png"))
	p.clf()
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(log_m200c, log_MF_c, c=data["redshift"][ok], s=5, marker='o',label="data", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r'log$_{10} n(>M)$') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	p.xlim((9.5,16))
	#p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-data.png"))
	p.clf()
	
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	sc1=p.scatter(logsigM1, log_f_c, c=data["redshift"][ok], s=5, marker='o',label="data", rasterized=True, vmin=zmin, vmax = zmax)
	sc1.set_edgecolor('face')
	cb = p.colorbar(shrink=0.8)
	cb.set_label("redshift")
	p.xlabel(r'$log(\sigma^{-1})$')
	p.ylabel(r'log$_{10} n(>M)$') # log$_{10}[ n(>M)]')
	gl = p.legend(loc=3,fontsize=10)
	gl.set_frame_on(False)
	#p.ylim((-8,1))
	#p.xlim((9.5,16))
	#p.xlim((-0.7,0.6))
	#p.ylim((-5.5,2))
	#p.yscale('log')
	p.grid()
	p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-data-xSigma.png"))
	p.clf()
	"""
plotData(qty = 'M200c', cos = "sat")
