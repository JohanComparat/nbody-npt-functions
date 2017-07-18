import StellarMass

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

import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()


smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
ll_dir = os.path.join(os.environ['GIT_NBODY_NPT'], 'data', 'stellar_mass_function')
path_ilbert13_SMF = os.path.join(ll_dir, "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(path_ilbert13_SMF, unpack=True)

smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
smf08 = lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )

mbins = n.arange(8,12.5,0.25)


import matplotlib.pyplot as p
out_dir = os.path.join(os.path.join(os.environ['MD10'],"results","stellar_mass_function", "images"))
	
# compare the stellar mass function measured to the Ilbert function
# take the AGN HGMF model

def plot_SMF_DC(snap_name, redshift):
	"""
	Plots the stellar mass functions and the corresponding duty cycle.
	"""
	# path for the output file
	out_duty_cycle = os.path.join(os.environ['MD10'],"duty_cycle", snap_name + "_duty_cycle.txt")
	# path for stellar mass function
	out_SMF_agn = os.path.join(os.environ['MD10'],"duty_cycle", snap_name + "_SMF.txt")
	out_SMF = os.path.join(os.environ['MD10'],"results", "stellar_mass_function", "data", "out_" + snap_name + "_SMF.txt")
	
	# path for duty cycle
	log_stellar_mass, duty_cycle = n.loadtxt(out_duty_cycle, unpack=True)
	dc = interp1d(log_stellar_mass, duty_cycle)
	p.figure(1, (6,6))
	p.plot(mbins, n.log10(smf01(10**mbins)), label='Ilbert 13, 0.2<z<0.5', ls='dashed')
	p.plot(mbins, n.log10(smf08(10**mbins)), label='Ilbert 13, 0.8<z<1.1', ls='dashed')
	print(out_SMF_agn)
	#logMs_low, logMs_up, counts, dN_dVdlogM = n.loadtxt(out_SMF_agn, unpack=True) 
	#p.plot((logMs_low+ logMs_up)/2., n.log10(dN_dVdlogM), label='MD10 AGN', lw=2)
	#print(out_SMF)
	
	logMs_low, logMs_up, counts, dN_dVdlogM_g = n.loadtxt(out_SMF, unpack=True) 
	p.plot((logMs_low+ logMs_up)/2., n.log10(dN_dVdlogM_g), label='MD10 GAL', lw=2)
	
	p.plot((logMs_low+ logMs_up)/2., n.log10(dc((logMs_low+ logMs_up)/2.)*dN_dVdlogM_g), label='MD10 AGN', lw=2)
	
	p.plot(mbins, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, redshift)) for logMs_i in mbins]) , label='Bo16')
	p.xlabel('stellar mass')
	p.ylabel('log Phi stellar mass')
	p.xlim((7., 12.2))
	p.ylim((-7,-1))
	p.title('z='+str(n.round(redshift,3)))
	p.grid()
	p.legend(loc=0, frameon=False)
	p.savefig(os.path.join(out_dir, "MD10_"+snap_name.zfill(5)+"_SMF.png"))
	p.clf()


# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for el in summ:
	print el
	plot_SMF_DC(el['snap_name'], el['redshift'])



#p.figure(1, (6,6))
#p.plot(logMS_DC_04, duty_cycle_04, label='MD 04')
#p.plot(logMS_DC_10, duty_cycle_10, label='MD 10')
#p.plot(logMS_DC_25, duty_cycle_25, label='MD 25')

#p.plot(logMS_DC_04_h, duty_cycle_04_h, label='MD h 04')
#p.plot(logMS_DC_10_h, duty_cycle_10_h, label='MD h 10')
#p.plot(logMS_DC_25_h, duty_cycle_25_h, label='MD h 25')

#p.axvline(7.2, c='k'  , ls='dashed')
#p.axvline(9.7, c='k' , ls='dashed')
#p.axvline(11.3, c='k', ls='dashed')
#p.xlabel('active fraction')
#p.ylabel('log stellar mass')
#p.xlim((6.5,12.2))
#p.yscale('log')
#p.ylim((0.005, .9))
#p.grid()
#p.legend(loc=0, frameon=False)
#p.savefig('/home/comparat/data/eRoMok/BO12_duty_cycle.png')
#p.clf()

