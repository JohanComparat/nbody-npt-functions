#import StellarMass

import numpy as n
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

bins = n.arange(6,13,0.1)

f_z = lambda z : n.piecewise(z, [z <= 1.1, z > 1.1], [ lambda z : (1.+z)**(5.82), lambda z : (1. + 1.1)**(5.82) * ((1.+z)/(1.+1.1))**(2.36)])

f_Mstar =lambda logM, z : (10**(logM - 10.99) )**(0.24) * n.e**( - 10**(logM - 10.99) )
	
def f_lambda_sar( logM, z, log_lambda_SAR ):
	lambda_SAR_var = 10**( log_lambda_SAR - 33.8 + 0.48 * (logM - 11.) )		
	g1z = 1.01 - 0.58 * (z - 1.1)
	return 1. / ( lambda_SAR_var**(g1z) + lambda_SAR_var**(3.72) )

psi_log = lambda log_lambda_SAR, logM, z : 10**(- 6.86) * f_lambda_sar( logM, z, log_lambda_SAR ) * f_Mstar(logM, z) * f_z(z)	

psi = lambda lambda_SAR, mass, z : psi_log(n.log10(lambda_SAR), n.log10(mass), z)

stellar_mass = 10**(10.8)
redshift=2.

integrand = lambda lambda_SAR, mass, z : psi( lambda_SAR, mass, z) / lambda_SAR / n.log(10)
fel = n.array([quad(integrand, 10**32, 10**33, args=(stellar_mass, redshift))[0], quad(integrand, 10**33, 10**34, args=(stellar_mass, redshift))[0], quad(integrand, 10**34, 10**35, args=(stellar_mass, redshift))[0], quad(integrand, 10**35, 10**36, args=(stellar_mass, redshift))[0]])
print(n.log10(fel), n.log10(n.sum(fel)))

import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()

import matplotlib.pyplot as p

# compare the stellar mass function measured to the Ilbert function
# take the AGN HGMF model

def tabulate_duty_cycle(snap_name, z, volume=1000.**3., bins = n.arange(6,13,0.1)):
	"""
	Tabulates the stellar mass functions in each snapshots and the corresponding duty cycle to be applied to recover the Bongiorno AGN function.
	"""
	# files to loop over
	fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+snap_name+'_SAM_Nb_?_Ms.fits')))
	fileList_snap.sort()
	print fileList_snap
	# path for the output file
	out_duty_cycle = os.path.join(os.environ['MD10'],"duty_cycle", "out_"+snap_name + "_duty_cycle.txt")
	# path for stellar mass function
	
	path_2_SMF = os.path.join(os.environ['MD10'], "duty_cycle", "out_"+snap_name+"_SMF.txt")
	# opens stellar mass function
	logMs_low, logMs_up, counts, dN_dVdlogM = n.loadtxt(path_2_SMF, unpack=True)
	# interpolates the duty cycle in the interesting region
	maxMS = n.max(logMs_up[(counts>1)])
	minMS = n.min(logMs_low[(counts>1)])
	x_SMF = (logMs_low+ logMs_up)/2.
	# interpolates the model of the AGN host galaxy mass function
	AGN_HGMF = n.array([xr.Phi_stellar_mass(logMs_i, z) for logMs_i in x_SMF])
	sel=(x_SMF>minMS)&(x_SMF<maxMS)
	duty_cycle = n.zeros_like(dN_dVdlogM)
	duty_cycle[sel] = AGN_HGMF[sel] / dN_dVdlogM[sel]
	duty_cycle[duty_cycle>1] = n.ones_like(duty_cycle[duty_cycle>1])
	print("HGMF", AGN_HGMF[sel])
	print("SMF", dN_dVdlogM[sel])
	print("DC",duty_cycle[sel],n.min(duty_cycle[sel]), n.max(duty_cycle[sel]))
	#duty_cycle = 10**AGN_HGMF(x_SMF[sel]) / dN_dVdlogM[sel]
	n.savetxt(out_duty_cycle, n.transpose([x_SMF, duty_cycle]), header = "log_stellar_mass duty_cycle")
    

# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for el in summ:
	print el
	tabulate_duty_cycle(el['snap_name'], el['redshift'])


