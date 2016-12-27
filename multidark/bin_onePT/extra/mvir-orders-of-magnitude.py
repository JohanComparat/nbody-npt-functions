from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import sys

import lib_functions_1pt as lib

import astropy.cosmology as co
cosmo = co.Planck13
import astropy.units as uu
from hmf import MassFunction

import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p


from scipy.interpolate import interp1d

from scipy.misc import derivative
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)

#Quantity studied
qty = "mvir"
# working directory
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits"))[1].data

NminCount = 100
Npmin = 300
nolim = [0,1e17]
limits_04 =  n.log10([Npmin*9.63 * 10**7, 5e12])
limits_10 =  n.log10([Npmin*1.51 * 10**9., 5e13])
limits_25 =  n.log10([Npmin*2.359 * 10**10., 5e14])
limits_40 =  n.log10([Npmin* 9.6 * 10**10. , 5e15])
MPART = n.array([9.63 * 10**7, 1.51 * 10**9, 2.359 * 10**10, 9.6 * 10**10])

zmin = -0.01
zmax = 2.3

#=======================
#=======================
cos = 'cen'
#=======================
#=======================

# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.mSelection(data, qty, limits_04, limits_10, limits_25,limits_40) 
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelCen)
# selection per box :
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc_new_rockS')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')


from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)

h0 = MassFunction(cosmo_model=cosmo, sigma_8=0.8229, z=0., delta_h=DeltaVir_bn98(0.), delta_wrt='mean')
h1 = MassFunction(cosmo_model=cosmo, sigma_8=0.8229, z=1., delta_h=DeltaVir_bn98(1.), delta_wrt='mean')
h2 = MassFunction(cosmo_model=cosmo, sigma_8=0.8229, z=2., delta_h=DeltaVir_bn98(2.), delta_wrt='mean')

hh = lambda zz : MassFunction(cosmo_model=cosmo, sigma_8=0.8229, z=zz, delta_h=DeltaVir_bn98(zz), delta_wrt='mean', Mmin=-5, Mmax=16.5)

delta_c_z = lambda z : h0.delta_c / hh(z).growth_factor

# spherical collaspe vs. redshift
zs = n.arange(0,2.5,0.1)
dc = n.empty_like(zs)
mc = n.empty_like(zs)
rc = n.empty_like(zs)
for ii, zz in enumerate(zs):
	h=hh(zz)
	dc[ii] = h0.delta_c / h.growth_factor
	s2m = interp1d(h.sigma, h.m)
	mc[ii] = s2m(dc[ii])
	rc[ii] = n.log10( (mc[ii]/(4*h.delta_h*n.pi*h.mean_density/3))**(1./3.)*1000.)
	print zz, dc[ii], rc[ii], n.log10(mc[ii])
	
fig, ax1 = p.subplots()
ax1.plot(zs, dc, 'k' )
ax1.set_xlabel('redshift')
ax1.set_ylabel(r'$\delta_{sc}(z)$')
ax1.set_yticks(dc[::5])
ax1.set_yticklabels(n.round(dc[::5],2))
ax1.grid()
ax2 = ax1.twinx()
ax2.plot(zs, dc, 'k' )
ax2.set_ylabel('log10(M collapsing)')
ax2.set_yticks(dc[::5])
ax2.set_yticklabels(n.round(n.log10(mc[::5]),2))
p.savefig(join(dir,"spherical-collapse.png"))
p.clf()

fig, ax1 = p.subplots()
ax1.plot(zs, rc, 'k' )
ax1.set_xlabel('redshift')
ax1.set_ylabel(r'log10($r_{vir}(z)$ [kpc])')
ax1.set_yticks(rc[::5])
ax1.set_yticklabels(n.round(rc[::5],2))
ax1.grid()
ax2 = ax1.twinx()
ax2.plot(zs, rc, 'k' )
ax2.set_ylabel('log10(M collapsing)')
ax2.set_yticks(rc[::5])
ax2.set_yticklabels(n.round(n.log10(mc[::5]),2))
p.savefig(join(dir,"spherical-collapse-r.png"))
p.clf()

n.savetxt(join("..","data","sc-relation.txt"), n.transpose([zs, dc, mc, rc]), header="z delta_c mvir_star rvir_star")

hrr = (cosmo.H(data['redshift']).value/100./cosmo.h)**3.
#h2.mean_density0

# x coordinates definition
logsig = n.log10(data['sigmaM'])#
lognu = n.log10(data['nuM'])
log_mvir = (data["log_"+qty+"_min"]+data["log_"+qty+"_max"])/2.
mvir = 10**log_mvir

# y coordinates
log_MF = n.log10( mvir * data["dNdVdlnM_"+cos]/ data["rhom"] )
log_MF_c = n.log10(  data["dNdVdlnM_"+cos+"_c"])

ff = mvir * data["dNdVdlnM_"+cos]/ data["rhom"]  / abs(data["dlnsigmaM1_o_dlnM"]) *hrr
ff_nu = data['nuM']*ff
log_f =  n.log10(ff)
log_f_nu =  n.log10(ff_nu)

log_f_c =  n.log10(mvir * data["dNdVdlnM_"+cos+"_c"]/ data["rhom"]  / abs(data["dlnsigmaM1_o_dlnM"]))

# NOW PLOTTING ALL THE DATA
#lib.plot_mvir_function_data(log_mvir[ok], logsigM1[ok], log_MF[ok], log_MF_c[ok], data['redshift'][ok], zmin, zmax, cos = cos)
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(lognu[ok], log_f[ok], c=data['redshift'][ok], s=5, marker='o',label="MD "+cos+" data", rasterized=True, vmin=zmin, vmax = zmax)
p.plot(n.log10(h0.nu), n.log10(h0.fsigma),'k--')
p.plot(n.log10(h1.nu), n.log10(h1.fsigma),'k--')
p.plot(n.log10(h2.nu), n.log10(h2.fsigma),'r--')
p.plot(n.log10(h0.nu), n.log10(lib.fnu_SMT(h0.nu, 0.322, 0.84, 0.3 )), 'm--', lw=2)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'$log(\nu)$')
p.ylabel(r'log$_{10} (M^2/\rho_m) dn(M)/dM$') 
 # log$_{10}[ n(>M)]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.ylim((-2.,0))
p.grid()
p.savefig(join(dir,"mvir-"+cos+"-differential-function-data-xNu.png"))
p.clf()


h0 = MassFunction(cosmo_model=cosmo, sigma_8=0.8229, z=0., delta_h=DeltaVir_bn98(0.), delta_wrt='mean', Mmin=7, Mmax=16.5)
h0a = MassFunction(cosmo_model=cosmo, sigma_8=0.8229*1.004**0.5, z=0., delta_h=DeltaVir_bn98(0.), delta_wrt='mean', Mmin=7, Mmax=16.5)
h0b = MassFunction(cosmo_model=cosmo, sigma_8=0.8229*1.01**0.5, z=0., delta_h=DeltaVir_bn98(0.), delta_wrt='mean', Mmin=7, Mmax=16.5)
h1 = MassFunction(cosmo_model=cosmo, sigma_8=0.8229, z=1., delta_h=DeltaVir_bn98(1.), delta_wrt='mean', Mmin=7, Mmax=16.5)
h2 = MassFunction(cosmo_model=cosmo, sigma_8=0.8229, z=2., delta_h=DeltaVir_bn98(2.), delta_wrt='mean', Mmin=7, Mmax=16.5)


fnu_SMT = lambda nu, Anorm, a, q : Anorm *a * (2./n.pi)**(0.5) *  ( 1 + (a**2*nu) **(-q) ) * a*nu**0.5 * n.e**( - a**2*nu / 2.)

model = interp1d(-n.log10(h0.sigma), n.log10(h0.fsigma))

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(-logsig[ok], log_f[ok]-model(-logsig[ok]), c=data['redshift'][ok], s=5, marker='o',label="MD "+cos+" data", rasterized=True, vmin=zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'$log(\sigma)$')
p.ylabel(r'log data - log model') 
 # log$_{10}[ n(>M)]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.xlim((-0.6,0.4))
p.ylim((-.5,.5))
p.grid()
p.savefig(join(dir,"mvir-"+cos+"-differential-function-discrepancy-xNu.png"))
p.clf()


p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(-logsig[ok], log_f[ok], c=data['redshift'][ok], s=5, marker='o',label="MD "+cos+" data", rasterized=True, vmin=zmin, vmax = zmax)
p.plot(-n.log10(h0.sigma), n.log10(h0.fsigma),'k--')
p.plot(-n.log10(h1.sigma), n.log10(h1.fsigma),'k--')
p.plot(-n.log10(h2.sigma), n.log10(h2.fsigma),'r--')
#p.plot(-logsig[ok], n.log10(fnu_SMT(10**(lognu[ok]**2.), 0.333, 0.794, 0.247 )), 'm+', lw=2)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'$log(\sigma)$')
p.ylabel(r'log$_{10} (M^2/\rho_m) dn(M)/dM$') 
 # log$_{10}[ n(>M)]')
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.xlim((-0.6,0.4))
p.ylim((-2.,0))
p.grid()
p.savefig(join(dir,"mvir-"+cos+"-differential-function-data-xNu.png"))
p.clf()


logmass = data['log_mvir'][ok]
mass = 10**data['log_mvir'][ok]
#dlnbin = dX / mass
m2sigma = interp1d(h0.M, h0.sigma )
sig = m2sigma( mass )
# m nu relation: nu = (delta_c / sigma_m)**2
m2nu = interp1d(h0.M, h0.nu )
nnu = m2nu( mass )
# jacobian
toderive = interp1d(n.log(h0.M), n.log(h0.sigma))
dlnsigmadlnm = derivative(toderive, n.log(mass) )


p.plot(n.log10(h0.m), h0.dndlnm*h0.m/h0.mean_density/h0.fsigma,'r+')
p.plot(data['log_mvir'][ok], abs(data['dlnsigmaMdlnM'][ok]),'b+')
p.savefig(join(dir,"sigma-m-derivative.png"))
p.clf()

p.plot(n.log10(h0.m), h0.sigma,'r--')
p.plot(n.log10(h0a.m), h0a.sigma,'k--')
p.plot(n.log10(h0b.m), h0b.sigma,'b--')
p.plot(data['log_mvir'][ok], sig, 'g+')

p.plot(data['log_mvir'][ok], data['sigmaM'][ok],'b+')
p.savefig(join(dir,"sigma-m-relation.png"))
p.clf()


lib.plot_mvir_function_data_perBox(log_mvir, log_MF, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos=cos)

sys.exit()
