import MultiDark as MD
import sys
import glob
import os
import numpy as n
import astropy.io.fits as fits
import astropy.units as uu

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
omega = lambda zz: cosmoMD.Om0*(1+zz)**3. / cosmoMD.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)
	
aexps = n.loadtxt(os.path.join(os.environ["MD10"], "hlists", "summary_all_aexp_hlists.list"), unpack=True)
	
out_name = os.path.join(os.environ["MD10"], 'hlists_MD_1.0Gpc.ascii')

aexps = n.array(aexps)
redshift = 1./aexps-1.
array_age = 10**9* cosmoMD.age(redshift).value
ids = n.argsort(redshift)

header = 'aexp redshift age_yr rho_crit delta_vir'

DATA = n.transpose([
	aexps[ids],  
	redshift[ids], 
	array_age[ids], 
	cosmoMD.critical_density(redshift[ids]).to(uu.solMass/(uu.kpc**3))*cosmoMD.h**2., 
	DeltaVir_bn98(redshift[ids]) ])


n.savetxt(out_name, DATA ,header = header)
