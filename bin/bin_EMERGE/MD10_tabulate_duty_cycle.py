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

# compare the stellar mass function measured to the Ilbert function
# take the AGN HGMF model

def tabulate_duty_cycle(snap_name, z, volume=1000.**3., bins = n.arange(6,13,0.1)):
	"""
	Tabulates the stellar mass functions in each snapshots and the corresponding duty cycle to be applied to recover the Bongiorno AGN function.
	"""
	# files to loop over
	fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+snap_name+'_SAM_Nb_?_Ms.fits')))
	fileList_snap.sort()
	print fileList_snap
	# path for the output file
	out_duty_cycle = os.path.join(os.environ['MD10'],"duty_cycle", snap_name + "_duty_cycle.txt")
	# path for stellar mass function
	out_SMF = os.path.join(os.environ['MD10'],"duty_cycle", snap_name + "_SMF.txt")
		
	Hall = n.zeros((len(fileList_snap), len(bins)-1))
	for ii, fileName in enumerate(fileList_snap):
		t0=time.time()
		# opens all relevant files
		hd = fits.open(fileName)[1].data
		print hd['stellar_mass_Mo13_mvir'][:10], n.min(hd['stellar_mass_Mo13_mvir'])
		Hall[ii], bb = n.histogram(hd['stellar_mass_Mo13_mvir'], bins=bins)

	counts = n.sum(Hall, axis=0) 
	dN_dVdlogM = counts*0.6777**3./(bins[1:]-bins[:-1])/volume/n.log(10)
	n.savetxt(out_SMF,  n.transpose([bins[:-1], bins[1:], counts, dN_dVdlogM]), header = "logMs_low logMs_up counts dN_dVdlogM")
	logMs_low, logMs_up = bins[:-1], bins[1:]
				
	# opens stellar mass function
	#logMs_low, logMs_up, counts, dN_dVdlogM = n.loadtxt(out_SMF, unpack=True)
	# interpolates the model of the AGN host galaxy mass function
	AGN_HGMF = interp1d(n.arange(5.5,13.5,0.01), n.array([n.log10(xr.Phi_stellar_mass(logMs_i, z)) for logMs_i in n.arange(5.5,13.5,0.01)]))
	# interpolates the duty cycle in the interesting region
	maxMS = n.max(logMs_up[(counts>1)])
	minMS = n.min(logMs_low[(counts>1)])
	x_SMF = (logMs_low+ logMs_up)/2.
	sel=(x_SMF>minMS)&(x_SMF<maxMS)
	duty_cycle = 10**AGN_HGMF(x_SMF[sel]) / dN_dVdlogM[sel]
	n.savetxt(out_duty_cycle, n.transpose([x_SMF[sel], duty_cycle]), header = "log_stellar_mass duty_cycle")
    

# open the output file_type
summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	

for el in summ:
	print el
	#if len(fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?_Ms.fits'))))>0 :
	tabulate_duty_cycle(el['snap_name'], el['redshift'])




sys.exit()


smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
ll_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'literature')
path_ilbert13_SMF = os.path.join(ll_dir, "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(os.path.join( ll_dir, "ilbert_2013_mass_function_params.txt"), unpack=True)

smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
smf08 = lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )

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

logMS_DC_10_h, duty_cycle_10_h = get_dc(glob.glob(os.path.join("..", "..", "data", "hlist*1.0*Gpc*_SMF.txt"))[0])
logMS_DC_04_h, duty_cycle_04_h = get_dc(glob.glob(os.path.join("..", "..", "data", "hlist*0.4Gpc*_SMF.txt"))[0])
logMS_DC_25_h, duty_cycle_25_h = get_dc(glob.glob(os.path.join("..", "..", "data", "hlist*2.5Gpc*_SMF.txt"))[0])

p.figure(1, (6,6))
p.plot(logMS_DC_04, duty_cycle_04, label='MD 04')
p.plot(logMS_DC_10, duty_cycle_10, label='MD 10')
p.plot(logMS_DC_25, duty_cycle_25, label='MD 25')

p.plot(logMS_DC_04_h, duty_cycle_04_h, label='MD h 04')
p.plot(logMS_DC_10_h, duty_cycle_10_h, label='MD h 10')
p.plot(logMS_DC_25_h, duty_cycle_25_h, label='MD h 25')

p.axvline(7.2, c='k'  , ls='dashed')
p.axvline(9.7, c='k' , ls='dashed')
p.axvline(11.3, c='k', ls='dashed')
p.xlabel('active fraction')
p.ylabel('log stellar mass')
p.xlim((6.5,12.2))
p.yscale('log')
p.ylim((0.005, .9))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig('/home/comparat/data/eRoMok/BO12_duty_cycle.png')
p.clf()

#minMS_MD04 = 7.2
#minMS_MD10 = 9.7
#minMS_MD25 = 11.3

n.savetxt(os.path.join("..", "..", "data", "duty_cycle_out_0.4Gpc_0.74230.txt"), n.transpose([logMS_DC_04,duty_cycle_04]), header="stellar_mass duty_cycle")

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
