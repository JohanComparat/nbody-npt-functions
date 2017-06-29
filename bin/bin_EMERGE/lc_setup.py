import glob
import sys
import astropy.io.fits as fits
import os
from os.path import join
# numerical modules
import numpy as n
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import curve_fit

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
cosmoDS = FlatLambdaCDM(H0=68.46*u.km/u.s/u.Mpc, Om0=0.298734, Ob0=0.046961)
from scipy.interpolate import interp1d
from scipy.integrate import quad

summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

L_box = 1000/0.6777

d_2_z = interp1d(summ['comoving_distance'],summ['redshift'])

shells = L_box * n.arange(1,4,1)

idz = n.searchsorted(summ['redshift'],d_2_z(shells))
print( d_2_z(shells),summ['redshift'][idz])

z_middle = (summ['redshift'][1:]+summ['redshift'][:-1])*0.5

z_mins = n.hstack((summ['redshift'][0], z_middle))
z_maxs = n.hstack((z_middle, summ['redshift'][-1]))
z_snap = summ['redshift']
dz = z_maxs - z_mins

# determine dz 
