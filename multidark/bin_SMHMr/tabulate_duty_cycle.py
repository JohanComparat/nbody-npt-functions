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


smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
ll_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'literature')
path_ilbert13_SMF = os.path.join(ll_dir, "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(os.path.join( ll_dir, "ilbert_2013_mass_function_params.txt"), unpack=True)

smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
smf08 = lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )


bins = n.arange(6,13,0.1)
xb = (bins[1:] + bins[:-1]) / 2.

# compare the stellar mass function measured to the Ilbert function
# take the AGN HGMF model




logMs = n.arange(5.5,13.5,0.01)
mbins = n.arange(8,12.5,0.25)

p.figure(1, (6,6))
    
p.plot(mbins, n.log10(smf01(10**mbins)), label='Ilbert 13, 0.2<z<0.5', ls='dashed')
p.plot(mbins, n.log10(smf08(10**mbins)), label='Ilbert 13, 0.8<z<1.1', ls='dashed')
logMs_low, logMs_up, counts, dN_dVdlogM = n.loadtxt(glob.glob(os.path.join("..", "..", "data", "out*1.0Gpc*_SMF.txt"))[0], unpack=True)
p.plot((logMs_low+ logMs_up)/2.-n.log10(0.6777), n.log10(dN_dVdlogM*0.6777**3./n.log(10)), label='MD10 out z=0.3', lw=2)
logMs_low, logMs_up, counts, dN_dVdlogM = n.loadtxt(glob.glob(os.path.join("..", "..", "data", "out*0.4Gpc*_SMF.txt"))[0], unpack=True)
p.plot((logMs_low+ logMs_up)/2.-n.log10(0.6777), n.log10(dN_dVdlogM*0.6777**3./n.log(10)), label='MD04 out z=0.3', lw=2)
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 0.3)) for logMs_i in logMs]) , label='z=0.3')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 0.9)) for logMs_i in logMs]) , label='z=0.9')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 1.2)) for logMs_i in logMs]) , label='z=1.2')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 1.8)) for logMs_i in logMs]) , label='z=1.8')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 2.5)) for logMs_i in logMs]) , label='z=2.5')
p.xlabel('stellar mass')
p.ylabel('log Phi stellar mass')
p.title('duty cycle')
p.xlim((7., 12.2))
p.ylim((-7,-1))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig('/home/comparat/data/eRoMok/BO12_MO13_duty_cycle.png')
p.clf()

AGN_HGMF = interp1d(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 1/0.74230-1.)) for logMs_i in logMs]))
  
def get_dc(path_to_SMF = glob.glob(os.path.join("..", "..", "data", "out*0.4Gpc*_SMF.txt"))[0]):
    logMs_low, logMs_up, counts, dN_dVdlogM = n.loadtxt(path_to_SMF, unpack=True)
    maxMS = n.max(logMs_up[(counts>1)])-n.log10(0.6777)
    minMS = n.min(logMs_low[(counts>1)])-n.log10(0.6777)
    x_SMF = (logMs_low+ logMs_up)/2.-n.log10(0.6777)
    sel=(x_SMF>minMS)&(x_SMF<maxMS)
    y_SMF = n.log10(dN_dVdlogM[sel]*0.6777**3./n.log(10))
    duty_cycle = 10**AGN_HGMF(x_SMF[sel]) / 10**y_SMF
    return x_SMF[sel], duty_cycle

logMS_DC_10, duty_cycle_10 = get_dc(glob.glob(os.path.join("..", "..", "data", "out*1.0*Gpc*_SMF.txt"))[0])
logMS_DC_04,duty_cycle_04 = get_dc(glob.glob(os.path.join("..", "..", "data", "out*0.4Gpc*_SMF.txt"))[0])
logMS_DC_25,duty_cycle_25 = get_dc(glob.glob(os.path.join("..", "..", "data", "out*2.5Gpc*_SMF.txt"))[0])

p.figure(1, (6,6))
p.plot(logMS_DC_04, duty_cycle_04)
p.plot(logMS_DC_10, duty_cycle_10)
p.plot(logMS_DC_25, duty_cycle_25)
p.axvline(7., 'k--')
p.axvline(9.5, 'k--')
p.axvline(11.2, 'k--')
p.xlabel('active fraction')
p.ylabel('log stellar mass')
p.xlim((6.5,12.2))
p.yscale('log')
p.ylim((0.005, .9))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig('/home/comparat/data/eRoMok/BO12_duty_cycle.png')
p.clf()


n.savetxt(os.path.join("..", "..", "data", "duty_cycle_0.4Gpc_0.74230.txt"), n.transpose([logMS_DC_04, duty_cycle_04]), header="stellar_mass duty_cycle")

sys.exit()

def measureSMF(env='MD04', volume=400.**3., file_type="out"):
	fileList = n.array(glob.glob(os.path.join(os.environ[env], "catalogs", file_type+"*.Ms.fits")))
	fileList.sort()
	print fileList
	Hall = n.zeros((len(fileList), len(bins)-1))
	for ii, fileN in enumerate(fileList):
		print fileN
		hd = fits.open(fileN)[1].data	
		Hall[ii], bb = n.histogram(hd['Mgal_mvir_Mo13'], bins=bins)
	
	counts =n.sum(Hall, axis=0)
	n.savetxt(os.path.join(os.environ[env], "results", os.path.basename(fileN)[:-5]+'_SMF.txt'),  n.transpose([bins[:-1], bins[1:], counts, counts/(bins[1:]-bins[:-1])/volume]), header = "logMs_low logMs_up counts dN_dVdlogM")


measureSMF(env='MD04', volume=400.**3.,  file_type="hlist")
measureSMF(env='MD04', volume=400.**3.,  file_type="out")

measureSMF(env='MD10', volume=1000.**3., file_type="hlist")
measureSMF(env='MD10', volume=1000.**3., file_type="out")

measureSMF(env='MD25', volume=2500.**3., file_type="hlist")
measureSMF(env='MD25', volume=2500.**3., file_type="out")
