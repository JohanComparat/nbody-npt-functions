import astropy.units as uu
import astropy.cosmology as co
aa = co.Planck13
import math as m
from scipy.integrate import quad

aah = co.FlatLambdaCDM(H0=100.0 *uu.km / (uu.Mpc *uu.s), Om0=0.307, Tcmb0=2.725 *uu.K, Neff=3.05, m_nu=[ 0.  ,  0. ,   0.06]*uu.eV, Ob0=0.0483)
rhom0 = aah.critical_density0.to(uu.solMass*uu.Mpc**-3).value
"""
#aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aah.critical_density0.to(uu.solMass*uu.Mpc**-3).value

vol = lambda zmin, zmax, area : (aah.comoving_volume(zmax)-aah.comoving_volume(zmin))*n.pi*area/129600.

volELG = n.log10(vol(0.6,1.,1400).value)
volQSO = n.log10(vol(0.9,2.2,7500).value)
volLya = n.log10(vol(2.1,3.5,7500).value)
n.log10(vol(0.9,1.6,24000).value)
n.log10(vol(0.9,2.2,24000).value)
n.log10(vol(0.1,0.4,8000).value)
n.log10(vol(0.75,1.4,2).value)
n.log10(vol(0.2,1.,6).value)
n.log10(vol(0.2,1.4,0.6).value)
"""
from scipy.interpolate import interp1d
import numpy as n
import matplotlib
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle
from os.path import join
from scipy.optimize import minimize

# limits at z0
Npmin = 1000
limits_04 = [Npmin*9.63 * 10**7, 5e12]
limits_10 = [Npmin*1.51 * 10**9., 5e13]
limits_25 = [Npmin*2.359 * 10**10., 5e14]
limits_40 = [Npmin* 9.6 * 10**10. , 5e15]
zmin = 0.
zmax = 5
NDecimal = 3

# defining directories :
dir = ".." #join("D:","data","MultiDark")
zl_04 = join(dir,"MD_0.4Gpc","output_SMD.list")
zl_10 = join(dir,"MD_1Gpc","output_MDPL.list")
zl_25 = join(dir,"MD_2.5Gpc","output_BigMD.list")
zl_40 = join(dir,"MD_4Gpc","output_HMD.list")

n0_04, a0_04, z0_04 = n.loadtxt(zl_04,unpack=True)
n0_10, z0_10, a0_10 = n.loadtxt(zl_10,unpack=True)
n0_25, a0_25, z0_25 = n.loadtxt(zl_25,unpack=True)
n0_40, z0_40, a0_40 = n.loadtxt(zl_40,unpack=True)

print n.max(z0_04[z0_04<3]), len(z0_04[z0_04<3])
print n.max(z0_10[z0_10<3]), len(z0_10[z0_10<3])
print n.max(z0_25[z0_25<3]), len(z0_25[z0_25<3])
print n.max(z0_40[z0_40<3]), len(z0_40[z0_40<3])

fig = p.figure(0,(13,6))
#fig.axes([0.15,0.3,0.8, 0.5])
ax1 = fig.add_subplot(111)
bins = n.arange(0.,1,0.01)
ax1.hist(n.log10(1+z0_04), bins=bins,label='SMD', histtype='step', lw=4)
ax1.hist(n.log10(1+z0_10), bins=bins,label='MDPL', histtype='step',lw=2)
ax1.hist(n.log10(1+z0_25), bins=bins,label='BigMD', histtype='step')
ax1.hist(n.log10(1+z0_40), bins=bins,label='HMD', histtype='step')
ax1.plot([n.log10(1+0.9),n.log10(1+2.2)],[1.7,1.7], label='QSO', lw=3)
ax1.plot([n.log10(1+0.6),n.log10(1+1.0)],[1.5,1.5], label='ELG', lw=3)
ax1.plot([n.log10(1+0.4),n.log10(1+0.8)],[1.3,1.3], label='LRG', lw=3)

ax1.grid()
xtik = n.hstack((n.arange(0,1.1,0.2), [1.5, 2, 3]))
ages = n.array([n.round(aa.age(xt).value) for xt in xtik ])
ax1.set_xticks(n.log10(1+xtik))
ax1.set_xticklabels(xtik)
ax1.set_xlabel('redshift')
ax1.set_ylabel('N snapshots / dlog(1+z)=0.01')
gl = p.legend(loc=0)
gl.set_frame_on(False)
ax2 = ax1.twiny()
#ax2.set_xlabel('age Gyr')
ax2.set_xticks(ax1.get_xticks())
ax2.set_xticklabels(ages)
ax2.set_xlabel('age (Gyr)')
#ax1.set_xlim((0,0.6))
#ax1.set_ylim((0,10))
p.savefig(join("..","presentationPlots","Nsnapshots.png"))
p.clf()

name, Lbox, Npart, Mp = n.loadtxt(join("..","existing-boxes.txt"), unpack=True, dtype=[('name', '<S8'), ('Lbox', '<i4'), ('Npart', '<i4'), ('Mpart', '<f4')])

nameS, volumeS, MhaloS = n.loadtxt(join("..","surveys.txt"), unpack=True, dtype=[('name', '<S12'),  ('volumeS', '<f4'), ('MhaloS', '<f4')])


logVol, massN100, massN10k, massN1M = n.loadtxt(join("..","M200c","volume-number.txt"), unpack=True)

NpH = 300

funPP = lambda volume, npart: 100*aa.h*rhom0* volume / npart 
vols = n.logspace(5,13,100)

fig = p.figure(1,(6,6))
#fig.axes([0.15,0.3,0.8, 0.5])
ax1 = fig.add_subplot(111)
ax1.plot(vols, funPP(vols, 4000**3.), 'm--', label=r'4000$^3$')
ax1.plot(vols, funPP(vols, 10000**3.), 'y--', label=r'10000$^3$')
ax1.plot(vols, funPP(vols, 40000**3.), 'c--', label=r'40000$^3$')
sel = (massN100>0)
ax1.plot(10**logVol[sel], 10**massN100[sel], label = 'HMF 100 halos') 
sel = (massN10k>0)
ax1.plot(10**logVol[sel], 10**massN10k[sel], label = 'HMF 10k halos') 
sel = (massN1M>0)
ax1.plot(10**logVol[sel], 10**massN1M[sel], label = 'HMF 1M halos') 

for ii, el in enumerate(name):
	#ax1.plot(Lbox[ii]**3., 300*10**Mp[ii], 'kx')
	ax1.plot(Lbox[ii]**3., NpH*10**Mp[ii], 'k+')
	#ax1.plot(Lbox[ii]**3., 300*10**Mp[ii], 'k_')
	#ax1.arrow(Lbox[ii]**3.,NpH*10**Mp[ii],0, NpH*10**Mp[ii]*2, fc='k', ec='k',head_width=Lbox[ii]**3.*0.9, head_length=NpH*10**Mp[ii]*1.1)
	#ax1.arrow(Lbox[ii]**3.,NpH*10**Mp[ii],-(Lbox[ii])**3./2., 0, fc='k', ec='k')
	ax1.annotate(el, xy=(Lbox[ii]**3., NpH*10**Mp[ii]),fontsize=10)#,rotation=45) #, xytext=(Lbox[ii]**3.*1.07, 1.07*NpH*10**Mp[ii])

for ii, el in enumerate(nameS):
	#print el
	p.plot(10**volumeS[ii], 10**MhaloS[ii], 'b^')
	ax1.annotate(el.replace('\\n','\n'), xy=(10**volumeS[ii], 10**MhaloS[ii]),color='b',fontsize=10)#,rotation=45)#, xytext=(10**volumeS[ii]/10, 100*10**MhaloS[ii]),color='b',fontsize=11,rotation=45)

ax1.set_xscale('log')
ax1.set_yscale('log')
p.grid()
ax1.set_xlabel(r'volume [Mpc$^{3}/h^{-3}$]')
ax1.set_ylabel(r'Halo mass resolved [$h^{-1}M_\odot$]')
ax1.set_ylim((3e8,1e15))
ax1.set_xlim((1e5, 1e13))
p.title(str(NpH)+' particles per halo')
gl = p.legend(loc=4,fontsize=10)
gl.set_frame_on(False)
p.savefig(join("..","presentationPlots","MassHalo-Volume.png"))
p.clf()

sys.exit()


fig = p.figure(1,(6,6))
#fig.axes([0.15,0.3,0.8, 0.5])
ax1 = fig.add_subplot(111)
for ii, el in enumerate(name):
	p.plot(Lbox[ii], 300*10**Mp[ii], 'bo')
	ax1.annotate(el, xy=(Lbox[ii], 300*10**Mp[ii]), xytext=(Lbox[ii]*1.07, 1.07*300*10**Mp[ii]))

ax1.set_xscale('log')
ax1.set_yscale('log')
p.grid()
ax1.set_xlabel('box length [Mpc/h]')
ax1.set_ylabel('Halo mass resolved [Msun/h]')
p.savefig(join("..","presentationPlots","Mass-Lbox.png"))
p.clf()



fig = p.figure(1,(6,6))
#fig.axes([0.15,0.3,0.8, 0.5])
ax1 = fig.add_subplot(111)
for ii, el in enumerate(name):
	p.plot(Lbox[ii]**3., Npart[ii]**3., 'bo')
	ax1.annotate(el, xy=(Lbox[ii]**3., Npart[ii]**3.), xytext=(Lbox[ii]**3.*1.07, 1.07*Npart[ii]**3.))

ax1.set_xscale('log')
ax1.set_yscale('log')
p.grid()
ax1.set_xlabel(r'box volume [Mpc$^{3}/h^{-3}$]')
ax1.set_ylabel('Number of particle')
p.savefig(join("..","presentationPlots","Npart-Volume.png"))
p.clf()


fig = p.figure(1,(6,6))
#fig.axes([0.15,0.3,0.8, 0.5])
ax1 = fig.add_subplot(111)
for ii, el in enumerate(name):
	p.plot(Lbox[ii]**3., 10**-Mp[ii], 'bo')
	ax1.annotate(el, xy=(Lbox[ii]**3., 10**-Mp[ii]), xytext=(Lbox[ii]**3.*1.07, 1.07*10**-Mp[ii]))

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim((1e-7,1e-12))
ax1.set_xlim((1e7, 1e13))
p.grid()
ax1.set_xlabel(r'box volume [Mpc$^{3}/h^{-3}$]')
ax1.set_ylabel('Particle mass [h/Msun]')
p.savefig(join("..","presentationPlots","mass-1-Volume.png"))
p.clf()

sys.exit()
