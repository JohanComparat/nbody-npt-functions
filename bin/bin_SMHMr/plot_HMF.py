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

mbins = n.arange(8,14.5,0.25)


import matplotlib.pyplot as p
out_dir = os.path.join(os.path.join(os.environ['MD10'],"results","mvir_mass_function", "images"))
	
# compare the stellar mass function measured to the Ilbert function
# take the AGN HGMF model

def plot_HMF_DC(snap_name, redshift):
	"""
	Plots the stellar mass functions and the corresponding duty cycle.
	"""
	# path for the output file
	# path for stellar mass function
	out_HMF = os.path.join(os.environ['MD10'],"results", "mvir_mass_function", "data", "out_" + snap_name + "_HMF.txt")
	# path to tracer SMFs
	out_file = lambda tracer_name : os.path.join(os.environ['MD10'],"results", "mvir_mass_function", "data", "out_"+snap_name+"_"+tracer_name+"_HMF.txt")

	p.figure(1, (6,6))
	
	logMs_low, logMs_up, counts, dN_dVdlogM_g = n.loadtxt(out_HMF, unpack=True) 
	ok = (dN_dVdlogM_g>0)
	p.plot((logMs_low[ok] + logMs_up[ok])/2., n.log10(dN_dVdlogM_g[ok]), label='MD10', lw=2, ls='dotted')
		
	def plot_tracer(tracer_name='4MOST_S5_BCG'):
		file_name = out_file(tracer_name )
		print file_name
		if os.path.isfile(file_name) :
			#print tracer_name
			logMs_low, logMs_up, counts, dN_dVdlogM_g = n.loadtxt(file_name  , unpack=True )
			ok = (dN_dVdlogM_g>0)
			p.plot((logMs_low[ok] + logMs_up[ok])/2., n.log10(dN_dVdlogM_g[ok]), label=tracer_name, ls='dashed', lw=0.75)

	plot_tracer("4MOST_S5_BCG" )
	plot_tracer("4MOST_S5_GAL" )
	plot_tracer("4MOST_S6_AGN" )
	plot_tracer("4MOST_S8_BG1" )
	plot_tracer("4MOST_S8_BG2" )
	plot_tracer("4MOST_S8_ELG" )
	plot_tracer("4MOST_S8_QSO" )

	
	p.ylabel(r'$\log_{10}(dN/dV/dlogM_{vir})$')
	p.xlabel(r'$\log_{10}(M_{vir})$')
	p.xlim((11., 14.5))
	p.ylim((-8,-2))
	p.title('z='+str(n.round(redshift,3)))
	p.grid()
	p.legend(loc=0, frameon=False)
	p.savefig(os.path.join(out_dir, "MD10_"+snap_name.zfill(5)+"_HMF_tracers.png"))
	p.clf()


# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for el in summ:
	print el
	plot_HMF_DC(el['snap_name'], el['redshift'])

os.system('cp $MD10/results/mvir_mass_function/images/*.png ~/wwwDir/eRoMok/mvir_mass_function/')


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

