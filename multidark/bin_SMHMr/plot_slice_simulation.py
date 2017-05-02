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


def get_slice(env='MD04', file_type="out", aexp='0.74230'):

	# parameters of the slice
	xmin, ymin, zmin = 0., 0., 0.,
	xmax, ymax, zmax = 60., 60., 30.
	
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "catalogs", file_type+"_"+aexp+"_*Xray.fits" )))
	fileList.sort()
	print fileList
	
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
		print "N points = ", len(xd[zone])
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

	print "N points total", len(x)

	def plot_slice(y, x, z, mass, LX_AGN, LX_cluster):
		agn_selection = (active)&(LX_AGN>42)
		cluster_selection = (LX_cluster>42)
		print "N AGN in plot=", len(y[agn_selection])
		print "N cluster in plot=", len(y[cluster_selection])
		p.figure(1, (9,9))
		p.plot(x[agn_selection], y[agn_selection], 'b+', label="AGN LX>42", rasterized=True)
		p.plot(x[cluster_selection], y[cluster_selection], 'ro', label="Cluster LX>42", rasterized=True)
		p.xlabel(r'$x [Mpc/h]$')
		p.ylabel(r'$y [Mpc/h]$')
		p.grid()
		p.yscale('log')
		p.xlim((xmin, xmax))
		p.ylim((ymin, ymax))
		p.legend(frameon=False)
		#p.title('Duty cycle 1%')
		p.savefig(os.path.join(os.environ[env], "results", file_type+"_"+aexp+'_xy_sim_slice.pdf'))
		p.clf()
		

	plot_slice(y, x, z, mass, LX_AGN, LX_cluster)

get_slice(env='MD04', file_type="out", aexp='0.74230')
os.system("cp $MD04/results/*.pdf ~/wwwDir/eRoMok/plots/MD_0.4Gpc/")
get_slice(env='MD10', file_type="out", aexp='0.74980')
os.system("cp $MD10/results/*.pdf ~/wwwDir/eRoMok/plots/MD_1.0Gpc/")
get_slice(env='MD25', file_type="out", aexp='0.75440')
os.system("cp $MD25/results/*.pdf ~/wwwDir/eRoMok/plots/MD_2.5Gpc/")
