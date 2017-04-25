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

bins = n.arange(10,15,0.1)
xb = (bins[1:] + bins[:-1]) / 2.

def create_plots(env='MD04', file_type="out"):
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
		cut420 = (hd['pid']=-1) & (hd['AGN']) & (hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'] > 42.)
		cut425 = (hd['pid']=-1) & (hd['AGN']) & (hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'] > 42.5)
		cut430 = (hd['pid']=-1) & (hd['AGN']) & (hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'] > 43.)
		cut435 = (hd['pid']=-1) & (hd['AGN']) & (hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'] > 43.5)
		Hall[ii], bb = n.histogram(hd['mvir'], bins=bins)
		H420[ii], bb = n.histogram(hd['mvir'][cut420], bins=bins)
		H425[ii], bb = n.histogram(hd['mvir'][cut425], bins=bins)
		H430[ii], bb = n.histogram(hd['mvir'][cut430], bins=bins)
		H435[ii], bb = n.histogram(hd['mvir'][cut435], bins=bins)
	
	all_halos_i = n.sum(Hall, axis=0)
	sel = (all_halos_i>0)
	all_halos = all_halos_i[sel].astype('float')
	                                                              
	p.figure(1, (6,6))
	p.fill_between(n.log10(AL12_mass), AL12_hod_up, AL12_hod_low, color='b', alpha=0.5, label='Allevato 12' )
	p.plot(xb, n.sum(H420, axis=0)[sel]/all_halos, label= 'L>42.0')
	p.plot(xb, n.sum(H425, axis=0)[sel]/all_halos, label= 'L>42.5')
	p.plot(xb, n.sum(H430, axis=0)[sel]/all_halos, label= 'L>43.0')
	p.plot(xb, n.sum(H435, axis=0)[sel]/all_halos, label= 'L>43.5')
	p.xlabel(r'$\log_{10} M_\odot$')
	p.ylabel(r'N AGN / N halo')
	p.grid()
	p.yscale('log')
	p.ylim((0.001,1.01))
	p.legend(frameon=False)
	p.title('Duty cycle 1%')
	p.savefig(os.path.join(os.environ[env], "results", os.path.basename(fileN)[:-5]+'_HOD_LX.pdf'))
	p.clf()
	
create_plots(env='MD04', file_type="hlist")
os.system("cp $MD04/results/*.pdf ~/wwwDir/eRoMok/plots/MD_0.4Gpc/")

create_plots(env='MD10', file_type="hlist")
os.system("cp $MD10/results/*.pdf ~/wwwDir/eRoMok/plots/MD_1.0Gpc/")

create_plots(env='MD25', file_type="hlist")
os.system("cp $MD25/results/*.pdf ~/wwwDir/eRoMok/plots/MD_2.5Gpc/")

sys.exit()

bins = n.arange(10,15,0.1)

def create_plots(env='MD04'):
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "catalogs", "*.Ms.fits")))
	print fileList
	for fileN in fileList:
		print fileN
		hd = fits.open(fileN)[1].data	
		cut420 = (hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'] > 42.)
		cut425 = (hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'] > 42.5)
		cut430 = (hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'] > 43.)
		cut435 = (hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'] > 43.5)
		p.figure(1, (6,6))
		p.hist(hd['mvir'][cut420], bins=bins, normed =True, label='L>42.0', histtype='step', cumulative = True)
		p.hist(hd['mvir'][cut425], bins=bins, normed =True, label='L>42.5', histtype='step', cumulative = True)
		p.hist(hd['mvir'][cut430], bins=bins, normed =True, label='L>43.0', histtype='step', cumulative = True)
		p.hist(hd['mvir'][cut435], bins=bins, normed =True, label='L>43.5', histtype='step', cumulative = True)
		p.xlabel(r'$\log_{10} M_\odot$')
		p.ylabel(r'normed cumulative distribution')
		p.grid()
		p.legend(frameon=False)
		p.savefig(os.path.join(os.environ[env], "results", os.path.basename(fileN)[:-5]+'_mvir_hist.pdf'))
        p.clf()
        
create_plots(env='MD04')
create_plots(env='MD10')
create_plots(env='MD25')

os.system("cp $MD04/results/*.pdf ~/wwwDir/eRoMok/plots/MD_0.4Gpc/")
os.system("cp $MD10/results/*.pdf ~/wwwDir/eRoMok/plots/MD_1.0Gpc/")
os.system("cp $MD25/results/*.pdf ~/wwwDir/eRoMok/plots/MD_2.5Gpc/")

