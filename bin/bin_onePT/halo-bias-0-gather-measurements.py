import numpy as n

from colossus.cosmology import cosmology
cosmology.setCosmology('planck15')

params = {'flat': True, 'H0': 67.77, 'Om0': 0.307115, 'Ob0': 0.048206, 'sigma8': 0.8228, 'ns': 0.9600, 'relspecies': False}
cosmology.setCosmology('myCosmo', params)


from colossus.lss import peaks
delta_c = peaks.collapseOverdensity(corrections = True, z = 0)

#delta_c = 1.686
nu = n.arange(0.2, 8.0, 0.01)
Masses = peaks.massFromPeakHeight(nu, 0.0)
sigma = ( delta_c / nu)

#p.semilogx(Masses,sigma,'k--')
#p.show()
#from colossus.lss import mass_function

import astropy.units as uu
from scipy.interpolate import interp1d
import matplotlib
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle
from os.path import join
from scipy.optimize import minimize
import scipy.fftpack as f
import time
#from hankel import SphericalHankelTransform
import os
#import lib_functions_1pt as lib
from hmf import MassFunction

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

import numpy as n
import astropy.io.fits as fits

version = 'v4'
qty = "clustering"
afactor=1.
redshift = 0.
top_dir = os.path.join(os.environ['OBS_REPO'], 'data_2017/MultiDark')
Rs, xiR  = n.loadtxt(join(top_dir,"PK_DM_CLASS","MD_z23_xi_2017.dat"), unpack = True)
linXi = interp1d(Rs,xiR)
#print Rs, xiR		
fileXi = n.array(glob.glob( join(top_dir, "MD_*Gpc*",  version, qty,"out_*30_xiR.pkl")))
fileXi.sort()

NminParticles = 1000.
rmin = 8.
rmax=20.

m2sigma = interp1d(Masses, sigma )

print( fileXi )

def get_hf(sigma_val=0.8228, boxRedshift=0., delta_wrt='mean'):
	"""
	Halo mass function model for the MultiDark simulation.
	"""
	#hf0 = MassFunction(cosmo_model=cosmo, sigma_8=sigma_val, z=boxRedshift)
	omega = lambda zz: cosmoMD.Om0*(1+zz)**3. / cosmoMD.efunc(zz)**2
	DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)
	print("DeltaVir", DeltaVir_bn98(boxRedshift), " at z",boxRedshift      )   
	hf1 = MassFunction(cosmo_model=cosmoMD, sigma_8=sigma_val, z=boxRedshift, delta_h=DeltaVir_bn98(boxRedshift), delta_wrt=delta_wrt, Mmin=8, Mmax=16.5)
	return hf1

hf = get_hf()

def get_basic_info(fileC, boxZN, delta_wrt='mean'):
	"""
	For a HMF measurement, this function returns all the basic information about the simulation used:
	- a HMF model with the corrected sigma8 value: hf,
	- size of the box: boxLength, boxLengthComoving
	- name of the simulation used later when writing results boxName
	- the redshiftof the measured HMF: boxRedshift, 7
	- the logmass of one particle: logmp,
	- correction to the mass measured due to force resolution: massCorrection
	"""
	if fileC.find('MD_0.4Gpc')>0:
		boxName='MD_0.4Gpc'
		nSN, aSN = n.loadtxt(join(os.environ['MD04_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		#hf = get_hf(0.8228*0.953**0.5, boxRedshift, delta_wrt=delta_wrt)
		# hf_ref = get_hf(0.8228, boxRedshift, delta_wrt=delta_wrt)
		logmp = n.log10(9.63 * 10**7/cosmoMD.h)
		boxLength = 400./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		massCorrection = 1. - 0.0002
		boxLengthComoving = 400.
		
	elif fileC.find('MD_1Gpc')>0 :
		boxName='MD_1Gpc'
		nSN, aSN = n.loadtxt(join(os.environ['MD10_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		#boxRedshift = 1./boxZN - 1.
		#hf = get_hf(0.8228*1.004**0.5, boxRedshift, delta_wrt=delta_wrt)
		# hf_ref = get_hf(0.8228, boxRedshift, delta_wrt=delta_wrt)
		logmp = n.log10(1.51 * 10**9/cosmoMD.h)
		boxLength = 1000./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		massCorrection = 1. - 0.0005
		boxLengthComoving = 1000.

	elif fileC.find('MD_2.5GpcNW')>0 :
		boxName='MD_2.5GpcNW'
		nSN, aSN = n.loadtxt(join(os.environ['MD25NW_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		#hf = get_hf(0.8228*1.01**0.5, boxRedshift, delta_wrt=delta_wrt)
		# hf_ref = get_hf(0.8228, boxRedshift, delta_wrt=delta_wrt)
		hz = 1.# hf.cosmoMD.H( boxRedshift ).value / 100.
		logmp = n.log10(2.359 * 10**10/cosmoMD.h )
		boxLength = 2500./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		massCorrection = 1. - 0.001
		boxLengthComoving = 2500.

	elif fileC.find('MD_4GpcNW')>0 :
		boxName='MD_4GpcNW'
		nSN, redshift40, aSN = n.loadtxt(join(os.environ['MD40NW_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'redshift', 'aSN'), 'formats': ('i4', 'f4', 'f4')})
		conversion = dict(n.transpose([ nSN, redshift40 ]))
		boxRedshift =  conversion[boxZN] 
		#hf = get_hf(0.8228*1.008**0.5, boxRedshift, delta_wrt=delta_wrt)
		# hf_ref = get_hf(0.8228, boxRedshift, delta_wrt=delta_wrt)
		hz = 1.# hf.cosmoMD.H( boxRedshift ).value / 100. 
		logmp = n.log10(9.6 * 10**10/cosmoMD.h )
		boxLength = 4000./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		boxLengthComoving = 4000.
		massCorrection = 1. - 0.003
	
	elif fileC.find('MD_2.5Gpc')>0 :
		boxName='MD_2.5Gpc'
		nSN, aSN, redshift25 = n.loadtxt(join(os.environ['MD25_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'aSN', 'redshift'), 'formats': ('i4', 'f4', 'f4')})
		conversion = dict(n.transpose([ nSN, redshift25 ]))
		boxRedshift =  conversion[boxZN] 
		#hf = get_hf(0.8228*1.01**0.5, boxRedshift, delta_wrt=delta_wrt)
		# hf_ref = get_hf(0.8228, boxRedshift, delta_wrt=delta_wrt)
		hz = 1.# hf.cosmoMD.H( boxRedshift ).value / 100.
		logmp = n.log10(2.359 * 10**10/cosmoMD.h)
		boxLength = 2500./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		boxLengthComoving = 2500.
		massCorrection = 1. - 0.001

	elif fileC.find('MD_4Gpc')>0 :
		boxName='MD_4Gpc'
		nSN, redshift40, aSN = n.loadtxt(join(os.environ['MD40_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'redshift', 'aSN'), 'formats': ('i4', 'f4', 'f4')})
		conversion = dict(n.transpose([ nSN, redshift40 ]))
		boxRedshift =  conversion[boxZN] 
		#hf = get_hf(0.8228*1.008**0.5, boxRedshift, delta_wrt=delta_wrt)
		# hf_ref = get_hf(0.8228, boxRedshift, delta_wrt=delta_wrt)
		hz = 1.# hf.cosmoMD.H( boxRedshift ).value / 100.
		logmp = n.log10(9.6 * 10**10 /cosmoMD.h )
		boxLength = 4000./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		boxLengthComoving = 4000.
		massCorrection = 1. - 0.003
	
	elif fileC.find('ds14_0')>0 :
		boxName='DS_8Gpc'
		boxRedshift = 0.
		#hf = get_hf_ds(0.8355, boxRedshift, delta_wrt=delta_wrt)
		# hf_ref = get_hf(0.8228, boxRedshift, delta_wrt=delta_wrt)
		hz = 1.# hf.cosmoMD.H( boxRedshift ).value / 100.
		logmp = n.log10(9.6 * 10**10 /cosmoDS.h )
		boxLength = 8000./cosmoMD.h/cosmoDS.efunc(boxRedshift)
		boxLengthComoving = 8000.
		massCorrection = 1. #- 0.003
	
	return boxLength, boxName, boxRedshift, logmp, boxLengthComoving, massCorrection 

def getBias(file):
	fileN = os.path.basename(file)
	boxN = file.split('/')[3]
	boxZN = fileN.split('_')[1]
	f=open(file,'r')
	bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbinsL, vbinsH = cPickle.load(f)
	f.close()
	
	rr = (bin_xi3D[1:] + bin_xi3D[:-1])/2.
	# print file
	m0=float(fileN.split('_')[3])
	m1=float(fileN.split('_')[4])
	
	m_low= 10**(m0)/cosmoMD.h
	m_high = 10**(m1)/cosmoMD.h
	m_mean = 10**((m0 + m1)/2.)/cosmoMD.h
	print('masses', m_low, m_mean, m_high)
	boxLength, boxName, boxRedshift, logmp, boxLengthComoving, massCorrection = get_basic_info(file, int(boxZN), delta_wrt='mean')
	print('logMP, nData', logmp, nD )
	if m_mean>NminParticles*10**logmp and nD>10000:
		sig_low = m2sigma( m_high)
		sig_high = m2sigma( m_low)
		sig_mean = m2sigma( m_mean)
		print('sigma', sig_low, sig_mean, sig_high, sig_low< sig_mean,sig_mean< sig_high )
		xi = DR*volume/(dV * pairCount) -1.
		ok = (rr>rmin)&(rr<rmax)
		scale = (n.min(rr[ok])+n.max(rr[ok]))/2.
		bias = n.mean(xi[ok]/linXi(rr[ok]))
		biasErr = n.std(xi[ok]/linXi(rr[ok]))
		#print [sig_low, sig_high, sig_mean, scale, bias, biasErr, volume, afactor]
		print('bias', bias )
		return [sig_low, sig_high, sig_mean, scale, bias, biasErr, volume, afactor, logmp]
		
	else:
		#print [-99.99, -99.99, -99.99, -99.99, -99.99, -99.99, -99.99, -99.99, -99.99]
		return [-99.99, -99.99, -99.99, -99.99, -99.99, -99.99, -99.99, -99.99, -99.99]

out=[]
for el in fileXi :
	got = getBias(el)
	if got == None :
		pass
	else:	
		out.append(got)

out= n.array(out)

sig_low, sig_high, sig_mean, scale, bias, biasErr, vol, aon, logmp = n.transpose(out)

kp = (sig_low!=-99.99)&(n.isnan(bias)==False)

tableOut = n.transpose([sig_low[kp], sig_high[kp], sig_mean[kp], scale[kp], bias[kp]**0.5, 0.5*biasErr[kp]*bias[kp]**(-0.5), vol[kp], aon[kp],logmp[kp]])

mvir_dir = join(os.environ['OBS_REPO'], 'data_2017/', 'mvir_dir')

n.savetxt(join(mvir_dir,  "halo-bias-measurement-summary.data"), tableOut, header=" sigma_low sigma_high sigma_mean scale bias biasErr volume a logmp")
