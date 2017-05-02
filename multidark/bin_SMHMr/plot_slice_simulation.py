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
	#agn_selection = (active)&(LX_AGN>42)
	#cluster_selection = (LX_cluster>42)
	
	zone = (selection_spatial) # & ( agn_selection | cluster_selection )

	return xd[zone], yd[zone], zd[zone], stellar_mass[zone], LX_AGN[zone], LX_cluster[zone], active[zone]

y, x, z, mass, LX_AGN, LX_cluster, active = get_plot_data(fileN[0])
for fileN in fileList[1:]:
	print fileN
	y_i, x_i, z_i, mass_i, LX_AGN_i, LX_cluster_i, active_i = get_plot_data(fileN)
	x = n.hstack((x, x_i))
	y = n.hstack((y, y_i))
	z = n.hstack((z, z_i))
	mass = n.hstack((mass, mass_i))
	LX_AGN = n.hstack((LX_AGN, LX_AGN_i))
	LX_cluster = n.hstack((LX_cluster, LX_cluster_i))
	active = n.hstack((active, active_i))


def plot_slice(y, x, z, mass, LX_AGN, LX_cluster):
	agn_selection = (active)&(LX_AGN>42)
	cluster_selection = (LX_cluster>42)
	p.figure(1, (6,6))
	p.plot(x[agn_selection], y[agn_selection], 'b+', label="AGN LX>42")
	p.plot(x[cluster_selection], y[cluster_selection], 'ro', label="Cluster LX>42")
	p.xlabel(r'$x /[Mpc/h]$')
	p.ylabel(r'$y /[Mpc/h]$')
	p.grid()
	p.yscale('log')
	p.ylim((0.000001,1.01))
	p.legend(frameon=False)
	#p.title('Duty cycle 1%')
	p.savefig(os.path.join(os.environ[env], "results", os.path.basename(fileN)[:-5]+'_sim_slice.pdf'))
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
