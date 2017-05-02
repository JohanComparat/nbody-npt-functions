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


#def get_slice(env='MD04', file_type="out", aexp='0.74230', out_dir = os.path.join("../../data/")):

env='MD04' 
file_type="out" 
aexp='0.74230' 
out_dir = os.path.join("../../data/")

# parameters of the slice
xmin, ymin, zmin = 0., 0., 0.,
xmax, ymax, zmax = 60., 60., 30.

# gets the file list to add the Xray luminosity
fileList = n.array(glob.glob(os.path.join(os.environ[env], "catalogs", file_type+"_"+aexp+"_*Xray.fits" )))

def get_plot_data(fileN):
	hd = fits.open(fileN)
	xd = hd[1].data['x']
	yd = hd[1].data['y']
	zd = hd[1].data['z']
	
	stellar_mass = hd[1].data['stellar_mass_Mo13_mvir']
	selection = hd[1].data['stellar_mass_reliable']

	LX_AGN = hd[1].data['lambda_sar_Bo16'] + hd[1].data['stellar_mass_Mo13_mvir']
	active = hd[1].data['activity']

	#hd[1].data['Mgas_cluster']
	#hd[1].data['kT_cluster']
	LX_cluster = hd[1].data['Lx_bol_cluster']
	#hd[1].data['Lx_ce_cluster']

	selection_spatial = (selection) & (xd > xmin) & (xd < xmax) & (yd > ymin) & (yd < ymax) & (zd > zmin) & (zd < zmax)
	agn_selection = (active)&(LX_AGN>42)
	cluster_selection = (LX_cluster>42)
	
	zone = (selection_spatial) # & ( agn_selection | cluster_selection )

	return xd[zone], yd[zone], zd[zone], stellar_mass[zone], LX_AGN[zone], LX_cluster[zone]

for fileN in fileList[1:]:
	print fileN
	hd = fits.open(fileN)[1].data	

	

def create_LF_plot(env='MD04', file_type="out"):
	
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "catalogs", file_type+"*.Ms.fits")))
	fileList.sort()
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
		defined = (hd['pid']==-1) & (hd['AGN']) 
		H420[ii], bb = n.histogram(hd['lambda_sar_Bo16'][defined]+hd['Mgal_mvir_Mo13'][defined], bins=bins)
	
	
	print env, volume_dict[env]                                                              
	p.figure(1, (6,6))
	p.fill_between(lglx[lf_sel], 10**phi_low[lf_sel], 10**phi_hi[lf_sel], color='b', alpha=0.5, label='Ueda 14' )
	
	all_halos_i = n.sum(Hall, axis=0)
	sel = (all_halos_i>0)
	all_halos = all_halos_i[sel].astype('float')
	p.plot(xb[sel], all_halos/volume_dict[env]/dlogX*n.log(10), label= 'simulation disctinct')
	
	all_halos_i = n.sum(H420, axis=0)
	sel = (all_halos_i>0)
	all_halos = all_halos_i[sel].astype('float')
	p.plot(xb[sel], all_halos/volume_dict[env]/dlogX*n.log(10), label= 'simulation AGN')
	
	p.xlabel(r'$\log_{10} (L_X/[erg/s])$')
	p.ylabel(r'luminosity function')
	p.grid()
	p.yscale('log')
	p.ylim((0.000001,1.01))
	p.legend(frameon=False)
	#p.title('Duty cycle 1%')
	p.savefig(os.path.join(os.environ[env], "results", os.path.basename(fileN)[:-5]+'_XLF.pdf'))
	p.clf()
	
create_LF_plot(env='MD04', file_type="hlist")
create_LF_plot(env='MD04', file_type="o")
os.system("cp $MD04/results/*.pdf ~/wwwDir/eRoMok/plots/MD_0.4Gpc/")

create_LF_plot(env='MD10', file_type="hlist")
create_LF_plot(env='MD10', file_type="out")
os.system("cp $MD10/results/*.pdf ~/wwwDir/eRoMok/plots/MD_1.0Gpc/")

create_LF_plot(env='MD25', file_type="hlist")
create_LF_plot(env='MD25', file_type="out")
os.system("cp $MD25/results/*.pdf ~/wwwDir/eRoMok/plots/MD_2.5Gpc/")
