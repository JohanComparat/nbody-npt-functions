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

from scipy.optimize import minimize
from scipy.optimize import curve_fit

from scipy.interpolate import interp1d
from scipy.misc import derivative

sigma = n.arange(0.05,10,0.05)
f = lambda sigma, A, a, b, c : A*( (sigma/b)**(-a) + 1 )*n.e**(-c/sigma**2.)
logFun = lambda logSigma, ps : n.log10( f(10.**logSigma, ps[0], ps[1], ps[2], ps[3]) )

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


errorLog = 0.03
NminCount = 1000
Npmin = 1000
limits_04 = [Npmin*9.63 * 10**7, 5e12]
limits_10 = [Npmin*1.51 * 10**9., 5e13]
limits_25 = [Npmin*2.359 * 10**10., 5e14]
limits_40 = [Npmin* 9.6 * 10**10. , 5e15]
MPART = n.array([9.63 * 10**7, 1.51 * 10**9, 2.359 * 10**10, 9.6 * 10**10])
names = n.array(["SMD", "MDPL", "BigMD", "HMD", "BigMDNW", "HMDNW"])

zmin = -0.01
zmax = 0.001
qty = 'M200c'
cos = "cen"
p0 = [0.224, 1.67, 1.80, 1.48]
p0 = [0.1, 1.2, 3.80, 1.3]
boundaries = [(0.01,0.1,0.1,0.1),(0.12,4,5,4)]
ps = p0
#def fitDataTinker08(qty = 'M200c', cos = "cen", zmin = zmin, zmax = zmax, p0 = [0.224, 1.67, 1.80, 1.48]):
"""
Plots the data to be used in the fits later in the analysis.
"""
# redshift selection
zSel = (data["redshift"]>zmin)&(data["redshift"]<zmax)
# mass selection
if  cos == "cen":
	mSel = ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0])) &(data["log_"+qty+"_max"]<n.log10(limits_04[1]))) | ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])) &(data["log_"+qty+"_max"]<n.log10(limits_10[1]))) |  ((data["boxLength"]==2500.)&(data["log_"+qty+"_min"]>n.log10(limits_25[0])) &(data["log_"+qty+"_max"]<n.log10(limits_25[1]))) |  ((data["boxLength"]==4000.)&(data["log_"+qty+"_min"]>n.log10(limits_40[0]))&(data["log_"+qty+"_max"]<n.log10(limits_40[1]))) 
if  cos == "sat":
	mSel =  ((data["boxLength"]==400.)&(data["log_"+qty+"_min"]>n.log10(limits_04[0])) &(data["log_"+qty+"_max"]<n.log10(limits_04[1]))) |  ((data["boxLength"]==1000.)&(data["log_"+qty+"_min"]>n.log10(limits_10[0])) &(data["log_"+qty+"_max"]<n.log10(limits_10[1])))
	
# minimum number counts selection
nSel = (data['dN_counts_'+cos]>NminCount)
# altogether
ok = (zSel) & (mSel) & (nSel)

# x coordinates
logsigM1 = n.log10(1./data['sigmaM'][ok])#
#print n.min(logsigM1), n.max(logsigM1)
log_m200c = (data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2.
m200c = 10**log_m200c
# mean density array normalization
rhom = cosmo.critical_density(data["redshift"][ok]).to(uu.solMass/(uu.Mpc)**3.)/(cosmo.H(data["redshift"][ok])/(100*uu.km/(uu.Mpc*uu.s)))**1.

x_data = n.log10(data['sigmaM'][ok])#
f_data = (m200c * data["dNdVdlnM_"+cos][ok]/ rhom.value  / abs(data["dlnsigmaM1_o_dlnM"][ok]) )
y_data = n.log10(f_data)
y_err = data['std90_pc_cen'][ok] / n.log(10.)

MD04=(ok) & (data["boxLength"]==400.)
MD10=(ok) & (data["boxLength"]==1000.)
MD25=(ok) & (data["boxLength"]==2500.)
MD40=(ok) & (data["boxLength"]==4000.)

qty = 'M200c'
cos = "cen"

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.plot(data["std90_pc_"+cos], data["dN_counts_"+cos]**(-0.5), 'ko', label='all', alpha=0.01)
p.plot(data["std90_pc_"+cos][MD04], data["dN_counts_"+cos][MD04]**(-0.5),marker='x',label="MD04",ls='')
p.plot(data["std90_pc_"+cos][MD10], data["dN_counts_"+cos][MD10]**(-0.5),marker='+',label="MD10",ls='')
p.plot(data["std90_pc_"+cos][MD25], data["dN_counts_"+cos][MD25]**(-0.5),marker='^',label="MD25",ls='')
p.plot(data["std90_pc_"+cos][MD40], data["dN_counts_"+cos][MD40]**(-0.5),marker='v',label="MD40",ls='')
xx = n.logspace(-4,0,20)
p.plot(xx, xx*3., ls='--', label='y=3x')
p.axhline(Npmin**-0.5, c='r', ls='--', label='min counts cut')#r'$1/\sqrt{10^3}$')
p.axhline((10**4.87)**-0.5, c='k', ls='--', label='min mass cut')#r'$1/\sqrt{10^{4.87}}$')
p.xlim((2e-4,4e-1))
p.ylim((2e-4,4e-1))
p.ylabel(r'$1/\sqrt{count} \; [\%]$')
p.xlabel(r'Jackknife  Resampling Error [%]')
p.yscale('log')
p.xscale('log')
gl = p.legend(loc=0,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join(dir,qty,"M200c-"+cos+"-jackknife-countsSqrt.png"))

p.show()

sys.exit()


# now the plot
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter((data["log_"+qty+"_min"][ok]+data["log_"+qty+"_max"][ok])/2., n.log10(data["std90_pc_"+cos+"_c"][ok]), c=data["redshift"][ok], s=5, marker='o',label="data", rasterized=True)
p.axhline(n.log10(0.05), label="5%")
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r'Error$[n(>M_{200c})]$ [%]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
#p.ylim((-8,1))
p.xlim((9.5,16))
#p.yscale('log')
p.grid()
p.savefig(join(dir,qty,"M200c-"+cos+"-cumulative-function-data-errPC.png"))
p.clf()
