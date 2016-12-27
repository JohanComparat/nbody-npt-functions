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

# defining directories :
dir = ".." #join("D:","data","MultiDark")
zList_file =  join(dir, "z-list-all-boxes.txt") 
z0 = n.loadtxt(zList_file,unpack=True)

# loads redshift 0 power spectrum
ks_lin_0, pks_lin_0 = n.loadtxt(join(dir, "Pk_DM_CLASS", "MD_z23_pk.dat"),unpack=True)
sel = (ks_lin_0<1e4)
ks_lin, pks_lin = ks_lin_0[sel], pks_lin_0[sel] *0.8228**2./0.7078
pk = interp1d(n.hstack((1e-20,ks_lin.min()/2., ks_lin, ks_lin.max()*2.,1e20)), n.hstack((0.,0.,pks_lin,0.,0.)))

dk = interp1d(pk.x, pk.x**3. * pk.y / (2*n.pi**2.) )
winTH = lambda k, R : 3*(m.sin(R*k) - R*k * m.cos(R*k)) / (R*k)**3.
integrand_sigma_2var = lambda lnk, R: dk(n.e**lnk) * (winTH(n.e**lnk,R))**2.

print pk.x, pk.y
# growth rate functionfrom Lahav et al. 91	
af = lambda z : 1 / (1+z)
g = lambda z : (5. * aa.Om(z) * af(z) / 2. ) / (aa.Om(z) - aa.Ode(z) + (1.+aa.Om(z)/2.)/(1.+aa.Ode(z)/70.))

#p.loglog(pk.x, pk.y,'b+')
#p.show()

def compute_M_sigma(zSnap, pk=pk, mmin=10**8.0, mmax=10**16.5):
	"""
	:param zSnap: redshift considered
	:param pk: power spectrum at redshift 0
	First computes critical density in (Msun/h) / (Mpc/h)^3
	Then computes radii and masses (M200c) to integrate the power spectrum over
	returns the radii, masses and sigma
	"""
	rho_cr = ( aa.critical_density(0.).to(uu.solMass/uu.megaparsec**3)*(aa.h**(-2)) ).value
	rmin =  n.log10( ( 3.*mmin / (200 * 4.*n.pi * rho_cr ) ) ** ( 1. / 3. ) )
	rmax = n.log10( ( 3.*mmax / (200 * 4.*n.pi * rho_cr ) ) ** ( 1. / 3. ) )
	Rs = n.logspace(rmin, rmax, 250) # Mpc/h
	print "redshift=", zSnap
	print "rho cr", rho_cr
	print "rmin rmax", rmin, rmax
	masses=200 * 4.*n.pi * rho_cr * (Rs)**3.0 / 3. # Msun/h
	print "masses", masses[0], masses[-1]
	sigma = n.empty(len(Rs))
	sigmaB = n.empty(len(Rs))
	for ii,R in enumerate(Rs):
		#integration in y = ln(k)
		integrand_sigma = lambda y: pk(n.e**y) * (m.sin(R*n.e**y) - R*n.e**y * m.cos(R*n.e**y))**2.0 / (n.e**y)**3.
		integr = quad(integrand_sigma,-20.7,20.7,limit=20001, epsrel=1.49e-08)  # correspond
		prefactor = 3. * g(zSnap) / ( 2**0.5 * n.pi * R ** 3.) 
		#integr = (3/(R**3.0) * (g(zSnap))**2. * quad(integrand_sigma,-9.0,9.,limit=20001, epsrel=1.49e-08) / (2* m.pi**2.0) # corresponds to sigma squared
		sigma[ii] = prefactor * integr[0]**0.5
		
		integrand_sigmaB = lambda lnk : integrand_sigma_2var(lnk, R)
		integrB = quad(integrand_sigmaB,-20.7,20.7,limit=20001, epsrel=1.49e-08)  # correspond
		sigmaB[ii] = integrB[0]**0.5

	print "sigma", n.min(sigma), n.max(sigma)
	return Rs,masses,sigma, sigmaB

zSnap=0.
Rs,masses,sigma,sigmaB = compute_M_sigma(zSnap)

p.plot(masses,sigma,'b+')
p.plot(masses,sigmaB,'r+')
#p.yscale('log')
p.xscale('log')
p.show()

sys.exit()

dk = interp1d(pk.x, pk.x**3. * pk.y / (2*n.pi**2.) )
winTH = lambda k, R : 3*(m.sin(R*k) - R*k * m.cos(R*k)) / (R*k)**3.
integrand_sigmaB = lambda lnk: dk(n.e**lnk)
integrB = quad(integrand_sigmaB,-20.7,20.7,limit=20001, epsrel=1.49e-08)  # correspond
quad(integrand_sigmaB,-20.7,0.,limit=20001, epsrel=1.49e-08)  # correspond

R=0.1
rho_cr = ( aa.critical_density(zSnap).to(uu.solMass/uu.megaparsec**3)*(aa.h**(-2)) ).value
mass=200 * 4.*n.pi * rho_cr * (8.0)**3.0 / 3.

R=8.
integrand_sigmaB = lambda lnk : integrand_sigma_2var(lnk, R) *0.8228**2./0.7078
integrB = quad(integrand_sigmaB,-20.7,20.7,limit=20001, epsrel=1.49e-08)  # correspond
sigmaB[ii] = integrB[0]**0.5


R=0.3
integrand_sigma = lambda y: pk(n.e**y) * 9.* (m.sin(R*n.e**y) - R*n.e**y * m.cos(R*n.e**y))**2.0 / (n.e**y *R )**3. / (2*n.pi**2.)
integr = quad(integrand_sigma,-20.7,20.7,limit=20001, epsrel=1.49e-08)   # 

prefactor = 3. * g(zSnap) / ( 2**0.5 * n.pi * R ** 3.) 	

sys.exit()

# loops over the redshifts of interest
for j, zSnap in enumerate(z0):
	Rs,masses,sigma = compute_M_sigma(zSnap)
	n.savetxt(join(dir,"Pk_DM_CLASS","MD_z"+str(zSnap)+"_Msigma.dat"),n.transpose([Rs,masses,sigma]), header=" R_hMpc M_hMsun sigma")
	

	
	