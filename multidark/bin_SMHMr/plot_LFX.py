import StellarMass
import XrayLuminosity

import numpy as n
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

print " set up box, and redshift "
#MD 1 hlist_0.74980_SAM_Nb_0.fits
#MD 25 hlist_0.75440_SAM_Nb_10.fits

#duty_cycle = 0.01

AL12_mass, AL12_hod_mean, AL12_hod_up, AL12_hod_low = n.loadtxt("/home/comparat/darksim/data/2PCF/allevato_2012_hod.txt", unpack=True)


z, lglx, phi_med, phi_low, phi_hi = n.loadtxt("/home/comparat/darksim/data/LXFunction/XLFAGN_STAN.RES". unpack=True)

#lf_sel = (z==0.25)
#0.25 38.5 -6.25795767806 -9.63144066207 -2.56060363285

dlogX = 0.2
bins = n.arange(38,45,dlogX)
xb = (bins[1:] + bins[:-1]) / 2.
volume_dict={'MD04': 400.**3., 'MD10': 1000.**3., 'MD25': 2500.**3.}

def create_LF_plot(env='MD04', file_type="out"):
	
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "catalogs", file_type+"*.Ms.fits")))
	print fileList
	Hall = n.zeros((len(fileList), len(bins)-1))
	H420 = n.zeros((len(fileList), len(bins)-1))
	H425 = n.zeros((len(fileList), len(bins)-1))
	H430 = n.zeros((len(fileList), len(bins)-1))
	H435 = n.zeros((len(fileList), len(bins)-1))
	for ii, fileN in enumerate(fileList):
		print fileN
		hd = fits.open(fileN)[1].data	
		defined = (hd['pid']==-1) #& (hd['AGN']) 
		Hall[ii], bb = n.histogram(hd['lambda_sar_Bo16'][defined]+hd['Mgal_mvir_Mo13'][defined], bins=bins)
	
	all_halos_i = n.sum(Hall, axis=0)
	sel = (all_halos_i>0)
	all_halos = all_halos_i[sel].astype('float')
	                                                              
	p.figure(1, (6,6))
	p.fill_between(lglx[lf_sel], phi_low[lf_sel], phi_hi[lf_sel], color='b', alpha=0.5, label='Ueda 14' )
	p.plot(xb, all_halos/volume_dict[env]/dlogX, label= 'simulation')
	p.xlabel(r'$\log_{10} (L_X/[erg/s])$')
	p.ylabel(r'luminosity function')
	p.grid()
	p.yscale('log')
	p.ylim((0.001,1.01))
	p.legend(frameon=False)
	p.title('Duty cycle 1%')
	p.savefig(os.path.join(os.environ[env], "results", os.path.basename(fileN)[:-5]+'_XLF.pdf'))
	p.clf()
	
create_plots(env='MD04', file_type="hlist")
os.system("cp $MD04/results/*.pdf ~/wwwDir/eRoMok/plots/MD_0.4Gpc/")

create_plots(env='MD10', file_type="hlist")
os.system("cp $MD10/results/*.pdf ~/wwwDir/eRoMok/plots/MD_1.0Gpc/")

create_plots(env='MD25', file_type="hlist")
os.system("cp $MD25/results/*.pdf ~/wwwDir/eRoMok/plots/MD_2.5Gpc/")
