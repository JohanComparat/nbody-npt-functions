from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import sys
import lib_functions_1pt as lib

import astropy.cosmology as co
cosmo = co.Planck13

#Quantity studied
qty = "vmax"
cos = 'cen'
# working directory
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits")
 )[1].data

NminCount = 100
limits_04 =  [125, 450] #max : 2e13#[100, 1000]
limits_10 =  [250, 800] #max : 2e14
limits_25 =  [600, 1100] #max : 5e14
limits_40 =  [900, 1400] #max : 2e14

zmin = 1. # -0.01
zmax = 2.3 # 1.0

# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.vSelection(data, qty, limits_04, limits_10, limits_25,limits_40) 
# minimum number counts selection
nSelCen = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelCen)
# selection per box :
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')


# x coordinates definition
vmax = data[qty]
log_vmax = n.log10(data[qty])

# y coordinates
norm = (100)**3. /(cosmo.H(data["redshift"]).value)**6.
log_VF = n.log10( norm * vmax**3. * data["dNdVdlnV_"+cos])
log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnV_"+cos+"_c"])
#print n.min(log_VF), n.max(log_VF)

# NOW PLOTTING ALL THE DATA
lib.plot_vmax_function_data(log_vmax[ok], log_VF[ok], log_VF_c[ok], data["redshift"][ok], zmin, zmax , cos=cos)

lib.plot_vmax_function_data_perBox(log_vmax, log_VF, log_VF_c, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos=cos)

"""
# PLOTTING THE ERROR PER BOX
lib.plot_vmax_function_data_error(log_vmax[ok & MD04], data['std90_pc_'+cos][ok & MD04], data["redshift"][ok & MD04], label='MD04', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data04-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD10], data['std90_pc_'+cos][ok & MD10], data["redshift"][ok & MD10], label='MD10', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data10-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25], data['std90_pc_'+cos][ok & MD25], data["redshift"][ok & MD25], label='MD25', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25NW], data['std90_pc_'+cos][ok & MD25NW], data["redshift"][ok & MD25NW], label='MD25NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25NW-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40], data['std90_pc_'+cos][ok & MD40], data["redshift"][ok & MD40], label='MD40', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40NW], data['std90_pc_'+cos][ok & MD40NW], data["redshift"][ok & MD40NW], label='MD40NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40NW-uncertainty.png")
"""
cos = 'sat'

# redshift selection
zSel = lib.zSelection( data, zmin, zmax )
# mass selection
mSel = lib.vSelection(data, qty, limits_04, limits_10, limits_25,limits_40) 
# minimum number counts selection
nSelSat = lib.nSelection(data, NminCount, cos )
# altogether
ok = (zSel) & (mSel) & (nSelSat)

# y coordinates
log_VF = n.log10( norm * vmax**3. * data["dNdVdlnV_"+cos])
log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnV_"+cos+"_c"])
#print n.min(log_VF), n.max(log_VF)

# NOW PLOTTING ALL THE DATA
lib.plot_vmax_function_data(log_vmax[ok], log_VF[ok], log_VF_c[ok], data["redshift"][ok], zmin , zmax , cos=cos)

lib.plot_vmax_function_data_perBox(log_vmax, log_VF, log_VF_c, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos=cos)

"""
# PLOTTING THE ERROR PER BOX
lib.plot_vmax_function_data_error(log_vmax[ok & MD04], data['std90_pc_'+cos][ok & MD04], data["redshift"][ok & MD04], label='MD04', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data04-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD10], data['std90_pc_'+cos][ok & MD10], data["redshift"][ok & MD10], label='MD10', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data10-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25], data['std90_pc_'+cos][ok & MD25], data["redshift"][ok & MD25], label='MD25', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD25NW], data['std90_pc_'+cos][ok & MD25NW], data["redshift"][ok & MD25NW], label='MD25NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data25NW-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40], data['std90_pc_'+cos][ok & MD40], data["redshift"][ok & MD40], label='MD40', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40-uncertainty.png")

lib.plot_vmax_function_data_error(log_vmax[ok & MD40NW], data['std90_pc_'+cos][ok & MD40NW], data["redshift"][ok & MD40NW], label='MD40NW', zmin = -0.01, zmax = 2.3, cos=cos, figName="vmax-"+cos+"-data40NW-uncertainty.png")

# ERROR PLOT: JK vs. POISSON
x = data["std90_pc_"+cos] 
y = data["dN_counts_"+cos]**(-0.5)
lib.plot_vmax_function_jackknife_poisson_error(x, y, MD04, MD10, MD25, MD25NW, MD40, MD40NW, cos = "cen")
"""