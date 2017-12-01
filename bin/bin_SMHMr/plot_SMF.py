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

import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()


smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
ll_dir = os.path.join(os.environ['GIT_NBODY_NPT'], 'data', 'stellar_mass_function')
path_ilbert13_SMF = os.path.join(ll_dir, "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(path_ilbert13_SMF, unpack=True)

smf_ilbert_fun = n.array([
lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
, lambda mass : smf_ilbert13( mass , 10**M_star[1], phi_1s[1]*10**(-3), alpha_1s[1], phi_2s[1]*10**(-3), alpha_2s[1] )
, lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )
, lambda mass : smf_ilbert13( mass , 10**M_star[3], phi_1s[3]*10**(-3), alpha_1s[3], phi_2s[3]*10**(-3), alpha_2s[3] )
, lambda mass : smf_ilbert13( mass , 10**M_star[4], phi_1s[4]*10**(-3), alpha_1s[4], phi_2s[4]*10**(-3), alpha_2s[4] )
, lambda mass : smf_ilbert13( mass , 10**M_star[5], phi_1s[5]*10**(-3), alpha_1s[5], phi_2s[5]*10**(-3), alpha_2s[5] )
, lambda mass : smf_ilbert13( mass , 10**M_star[6], phi_1s[6]*10**(-3), alpha_1s[6], phi_2s[6]*10**(-3), alpha_2s[6] )
, lambda mass : smf_ilbert13( mass , 10**M_star[7], phi_1s[7]*10**(-3), alpha_1s[7], phi_2s[7]*10**(-3), alpha_2s[7] )
])

mbins = n.arange(8,12.5,0.25)

smf_ilbert_zmin = n.array([ 
0.2
, 0.5
, 0.8
, 1.1
, 1.5
, 2.0
, 2.5
, 3.0 ])

smf_ilbert_zmax = n.array([ 
0.5
, 0.8
, 1.1
, 1.5
, 2.0
, 2.5
, 3.0
, 4.0 ])

smf_ilbert_name = n.array([ "Il13 "+str(zmin)+"<z<"+str(zmax) for zmin, zmax in zip(smf_ilbert_zmin,smf_ilbert_zmax) ])



import matplotlib.pyplot as p
out_dir = os.path.join(os.path.join(os.environ['MD10'],"results","stellar_mass_function", "images"))
	
# compare the stellar mass function measured to the Ilbert function
# take the AGN HGMF model

def plot_SMF_DC(snap_name, redshift):
	"""
	Plots the stellar mass functions and the corresponding duty cycle.
	"""
	# path for the output file
	out_duty_cycle = os.path.join(os.environ['MD10'],"duty_cycle", "out_"+snap_name + "_duty_cycle.txt")
	# path for stellar mass function
	out_SMF = os.path.join(os.environ['MD10'],"duty_cycle", "out_" + snap_name + "_SMF.txt")
	# path to tracer SMFs
	out_file = lambda tracer_name : os.path.join(os.environ['MD10'], "duty_cycle", "out_"+snap_name+"_"+tracer_name+"_SMF.txt")

	p.figure(1, (6,6))
	for fun, name in zip(smf_ilbert_fun, smf_ilbert_name):
		p.plot(mbins, n.log10(fun(10**mbins)), label=name, ls='dashed', lw=0.5)
	
	#p.plot(mbins, n.log10(smf08(10**mbins)), label='Ilbert 13, 0.8<z<1.1', ls='dashed')
	#print(out_SMF_agn)
	#logMs_low, logMs_up, counts, dN_dVdlogM = n.loadtxt(out_SMF_agn, unpack=True) 
	#p.plot((logMs_low+ logMs_up)/2., n.log10(dN_dVdlogM), label='MD10 AGN', lw=2)
	#print(out_SMF)
	
	# galaxies
	logMs_low, logMs_up, counts, dN_dVdlogM_g = n.loadtxt(out_SMF, unpack=True) 
	print logMs_low, dN_dVdlogM_g 
	ok = (counts>2)
	print "SMF", n.min(logMs_low[ok]), n.max(logMs_up[ok]) 
	p.plot((logMs_low[ok] + logMs_up[ok])/2., n.log10(dN_dVdlogM_g[ok]), label='MD10 GAL')#, lw=2)
	
	# path for duty cycle
	log_stellar_mass, duty_cycle = n.loadtxt(out_duty_cycle, unpack=True)
	print "DC",n.min(log_stellar_mass), n.max(log_stellar_mass)
	dc = interp1d(log_stellar_mass, duty_cycle)
	p.plot((logMs_low[ok]+ logMs_up[ok])/2., n.log10(dc((logMs_low[ok] + logMs_up[ok])/2.)*dN_dVdlogM_g[ok]), label='MD10 AGN')#, lw=2)
	
	#p.plot(mbins, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, redshift)) for logMs_i in mbins]) , label='Bo16', ls='dashed')
	
	def plot_tracer(tracer_name='4MOST_S5_BCG'):
		file_name = out_file(tracer_name )
		print file_name
		if os.path.isfile(file_name) :
			#print tracer_name
			logMs_low, logMs_up, counts, dN_dVdlogM_g = n.loadtxt(file_name  , unpack=True )
			ok = (dN_dVdlogM_g>0)
			p.plot((logMs_low[ok] + logMs_up[ok])/2., n.log10(dN_dVdlogM_g[ok]), label=tracer_name)#, ls='dashed', lw=0.75)

	#plot_tracer("4MOST_S5_BCG" )
	#plot_tracer("4MOST_S5_GAL" )
	#plot_tracer("4MOST_S6_AGN" )
	#plot_tracer("4MOST_S8_BG1" )
	#plot_tracer("4MOST_S8_BG2" )
	#plot_tracer("4MOST_S8_ELG" )
	#plot_tracer("4MOST_S8_QSO" )

	
	p.xlabel('stellar mass')
	p.ylabel('log Phi stellar mass')
	p.xlim((9., 12.2))
	p.ylim((-8,-2))
	p.title('z='+str(n.round(redshift,3)))
	p.grid()
	p.legend(loc=0, frameon=False)
	p.savefig(os.path.join(out_dir, "MD10_"+snap_name.zfill(5)+"_SMF_tracers.png"))
	p.clf()


# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for el in summ:
	print el
	plot_SMF_DC(el['snap_name'], el['redshift'])

#os.system('cp $MD10/results/stellar_mass_function/images/*.png ~/wwwDir/eRoMok/stellar_mass_function/')


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

