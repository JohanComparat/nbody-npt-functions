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

import matplotlib.pyplot as p

# stellar mass function measured to the Ilbert function
smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
ll_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'literature')
path_ilbert13_SMF = os.path.join(ll_dir, "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(os.path.join( ll_dir, "ilbert_2013_mass_function_params.txt"), unpack=True)

smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
smf08 = lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )

# the AGN HGMF model
def plot_duty_cycle(env='MD04', volume=400.**3.,  file_type="hlist", aexp='0.74230', out_dir = os.path.join("../../data/")):
    # path for the output file
    path_to_duty_cycle = os.path.join(out_dir, env+"_"+file_type+"_"+aexp+"_duty_cycle.txt")
    log_stellar_mass, duty_cycle = n.loadtxt(path_to_duty_cycle, unpack="True")
    p.plot(log_stellar_mass, duty_cycle, label=env+file_type+' a='+aexp)
    
def plot_SMF(env='MD04', volume=400.**3.,  file_type="hlist", aexp='0.74230', out_dir = os.path.join("../../data/")):
    # path for the output file
    path_to_SMF = os.path.join(out_dir, env+"_"+file_type+"_"+aexp+"_SMF.txt")
    logMs_low, logMs_up, counts, dN_dVdlogM = n.loadtxt(path_to_SMF, unpack=True)
    p.plot((logMs_low + logMs_up)/2., n.log10(dN_dVdlogM), label=env+' a='+aexp)

logMs = n.arange(7,12.5,0.25)

p.figure(1, (6,6))
# Ilbert 2013   
p.plot(logMs, n.log10(smf01(10**logMs)), label='Il13, 0.2<z<0.5', ls='dashed')
p.plot(logMs, n.log10(smf08(10**logMs)), label='Il13, 0.8<z<1.1', ls='dashed')
# MultiDark
#plot_SMF(env='MD04', volume=400.**3.,  file_type="hlist", aexp='0.74230', out_dir = os.path.join("../../data/"))
plot_SMF(env='MD04', volume=400.**3.,  file_type="out"  , aexp='0.74230', out_dir = os.path.join("../../data/"))
#plot_SMF(env='MD10', volume=1000.**3., file_type="hlist", aexp='0.74980', out_dir = os.path.join("../../data/"))
plot_SMF(env='MD10', volume=1000.**3., file_type="out"  , aexp='0.74980', out_dir = os.path.join("../../data/"))
#plot_SMF(env='MD25', volume=2500.**3., file_type="hlist", aexp='0.75440', out_dir = os.path.join("../../data/"))
plot_SMF(env='MD25', volume=2500.**3., file_type="out"  , aexp='0.75440', out_dir = os.path.join("../../data/"))
# Bongiorno 2012
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 0.3)) for logMs_i in logMs]) , label='BO13 z=0.3', ls='dotted')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 0.9)) for logMs_i in logMs]) , label='BO13 z=0.9', ls='dotted')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 1.8)) for logMs_i in logMs]) , label='BO13 z=2.0', ls='dotted')
p.xlabel(r'$\log(M_*/M_\odot)$')
p.ylabel(r'$\log(\Phi(M_*) / [Mpc^3 dex])$')
p.xlim((7., 12.5))
p.ylim((-7,-1))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig('/home/comparat/data/eRoMok/BO12_MO13_Il13_SMF.png')
p.clf()


p.figure(1, (6,6))
#plot_duty_cycle(env='MD04', volume=400.**3.,  file_type="hlist", aexp='0.74230', out_dir = os.path.join("../../data/"))
plot_duty_cycle(env='MD04', volume=400.**3.,  file_type="out"  , aexp='0.74230', out_dir = os.path.join("../../data/"))
#plot_duty_cycle(env='MD10', volume=1000.**3., file_type="hlist", aexp='0.74980', out_dir = os.path.join("../../data/"))
plot_duty_cycle(env='MD10', volume=1000.**3., file_type="out"  , aexp='0.74980', out_dir = os.path.join("../../data/"))
#plot_duty_cycle(env='MD25', volume=2500.**3., file_type="hlist", aexp='0.75440', out_dir = os.path.join("../../data/"))
plot_duty_cycle(env='MD25', volume=2500.**3., file_type="out"  , aexp='0.75440', out_dir = os.path.join("../../data/"))
p.axvline(7.2, c='k'  , ls='dashed')
p.axvline(9.7, c='k' , ls='dashed')
p.axvline(11.3, c='k', ls='dashed')
p.ylabel('active fraction [%]')
p.xlabel(r'$\log(M_*/M_\odot)$')
p.xlim((6.5,12.2))
p.yscale('log')
p.ylim((0.005, .9))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig('/home/comparat/data/eRoMok/BO12_duty_cycle.png')
p.clf()


