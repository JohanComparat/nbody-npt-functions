import astropy.units as uu
import astropy.cosmology as co
aa = co.Planck13
import math as m
from scipy.integrate import quad

#aah = co.FlatLambdaCDM(H0=100.0 *uu.km / (uu.Mpc *uu.s), Om0=0.307, Tcmb0=2.725 *uu.K, Neff=3.05, m_nu=[ 0.  ,  0. ,   0.06]*uu.eV, Ob0=0.0483)
#rhom = aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aah.critical_density0.to(uu.solMass*uu.Mpc**-3).value

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
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")

zList_file =  join(dir_25, "redshift_list.txt") 
n0,z0,a0 = n.loadtxt(zList_file,unpack=True)
n1 = n0.astype(int)+1
n2 = n1.astype(str)
snap2Z = {n2[i].zfill(3): z0[i] for i in range(len(n2))}

pklist = n.array(glob.glob(join(dir_25,"Pk_DM","PkDM_BigMD_Planck_snap_*.dat")))
"""
j=79
pkfile = pklist[j]
snapNum=pkfile.split('_')[-1][:-4] #"080"
zSnap = snap2Z[snapNum]
ks, pks = n.loadtxt(pkfile,unpack=True,usecols =(1,3))
p.loglog(ks, pks,label=str(n.round(zSnap,3)))
ks_lin_0, pks_lin_0 = n.loadtxt(join(dir,"Pk_DM_CLASS","multidark00_z"+n2[j]+"_pk.dat"),unpack=True)
sel = (ks_lin_0<1e4)
ks_lin, pks_lin = ks_lin_0[sel], pks_lin_0[sel]
p.loglog(ks_lin, pks_lin,label="CLASS")

p.ylim(1,3e4)
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.show()
"""
for j in n.arange(80):
	pkfile = pklist[j]
	snapNum=pkfile.split('_')[-1][:-4] #"080"
	zSnap = snap2Z[snapNum]
	ks_lin_0, pks_lin_0 = n.loadtxt(join(dir,"Pk_DM_CLASS","multidark00_z"+n2[j]+"_pk.dat"),unpack=True)
	sel = (ks_lin_0<1e4)
	ks_lin, pks_lin = ks_lin_0[sel], pks_lin_0[sel]
	pk = interp1d(n.hstack((1e-20,ks_lin.min()/2., ks_lin, ks_lin.max()*2.,1e20)), n.hstack((0.,0.,pks_lin,0.,0.)))
	rho_bar = (aa.critical_density(zSnap).to(uu.solMass/uu.megaparsec**3)*(aa.h**-2)).value
	Rs = n.logspace(-0.7,1.3,100)
	masses=4.*n.pi * rho_bar * (Rs)**3.0 / 3.
	sigma = n.empty(len(Rs))
	for ii,R in enumerate(Rs):
		integrand_sigma = lambda lk: pk((10**lk)) * ( (3/((10**lk)*R)**3.0) * (m.sin((10**lk)*R)-(10**lk)*R*m.cos((10**lk)*R)) )**2.0 * (10**lk) **3.0 *n.log(10)/ (2* m.pi**2.0)
		integr=quad(integrand_sigma,-9.0,9.,limit=20001, epsrel=1.49e-08) # corresponds to sigma squared
		sigma[ii] = integr[0]**0.5

	n.savetxt(join(dir,"Pk_DM_CLASS","multidark00_z"+n2[j]+"_Msigma.dat"),n.transpose([Rs,masses,sigma]), header=" R M sigma")
	print n2[j]
	
sys.exit()

for j in n.arange(79)[9:][::10]:
	pkfile = pklist[j]
	snapNum=pkfile.split('_')[-1][:-4] #"080"
	zSnap = snap2Z[snapNum]
	ks,pks = n.loadtxt(pkfile,unpack=True,usecols =(1,4))
	p.loglog(ks, pks,label=str(n.round(zSnap,3)))
	ks_class, pks_class = n.loadtxt(join(dir,"Pk_DM_CLASS","multidark00_z"+n2[j]+"_pk.dat"),unpack=True)
	p.loglog(ks_class, pks_class, ls='--')#,label="CLASS")
	print pklist[j], "multidark00_z"+n2[j]+"_pk.dat", zSnap
	
p.xlim((1e-3,1))
p.ylim(1,3e4)
gl = p.legend(loc=3,fontsize=10)
gl.set_frame_on(False)
p.show()
