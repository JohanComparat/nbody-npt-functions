import astropy.io.fits as fits
import matplotlib.pyplot as p
import numpy as n
from os.path import join
import os
import sys 

from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import sys

import lib_functions_1pt as lib
from lib_functions_1pt import *


from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
sigma_val=0.8229
delta_c = 1.686
from scipy.interpolate import interp1d
from scipy.integrate import quad
import numpy as n
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

boxRedshift = 0.
version='v3'

omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)
hf = MassFunction(cosmo_model=cosmo, sigma_8=sigma_val, z=boxRedshift, delta_h=DeltaVir_bn98(boxRedshift), delta_wrt='mean', Mmin=7, Mmax=16.5)

f_BH = lambda sigma, A, a, p, q: A* (2./n.pi)**(0.5) * ( 1 + (sigma**2./(a**delta_c*2.))**(p) )*(delta_c*a**0.5/sigma)**(q)*n.e**(-a*delta_c**2./(2.*sigma**2.))

X = n.arange(-0.6, 0.5, 0.01) #n.log10(1./sigma)
sigma = 10**-X 

hz = cosmo.H( boxRedshift ).value / 100.
# m sigma relation using the sigma8 corrected power spectrum
m2sigma = interp1d(hf.M, hf.sigma )
# m nu relation: nu = (delta_c / sigma_m)**2
m2nu = interp1d(hf.M, hf.nu )
# jacobian
toderive = interp1d(n.log(hf.M), n.log(hf.sigma))
mass=hf.M[100:-100]
dlnsigmadlnm = derivative(toderive, n.log(mass) )
rhom_units = cosmo.Om(boxRedshift)*cosmo.critical_density(boxRedshift).to(u.solMass/(u.Mpc)**3.)#/(cosmo.h)**2.
# in units (Msun/h) / (Mpc/h)**3
rhom = rhom_units.value # hf.mean_density#/(hz)**2. 

ftC16 = f_BH(hf.sigma[100:-100], 0.279, 0.908, 0.671, 1.737)

MF_MD = interp1d(mass, ftC16*rhom*abs(dlnsigmadlnm)/mass)

NpartMin = 50.
p_init = (-1.85, 7., -2.3, 4.)

mp04 = n.log10(NpartMin*9.63 * 10**7)
mp10 = n.log10(NpartMin*1.51 * 10**9)
mp25 = n.log10(NpartMin*2.359 * 10**10)
mp40 = n.log10(NpartMin*9.6 * 10**10. )

flist_1 = n.hstack(( n.array(glob.glob(join(os.environ['MD04_DIR'],version, "subhalos", "out_*_subhalos_inDistinct.fits"))), n.array(glob.glob(join(os.environ['MD10_DIR'],version, "subhalos", "out_*_subhalos_inDistinct.fits"))), n.array(glob.glob(join(os.environ['MD25_DIR'],version, "subhalos", "out_*_subhalos_inDistinct.fits"))), n.array(glob.glob(join(os.environ['MD25NW_DIR'],version, "subhalos", "out_*_subhalos_inDistinct.fits"))), n.array(glob.glob(join(os.environ['MD40_DIR'],version, "subhalos", "out_*_subhalos_inDistinct.fits"))), n.array(glob.glob(join(os.environ['MD40NW_DIR'],version, "subhalos", "out_*_subhalos_inDistinct.fits"))) ))
flist_2 = n.hstack(( n.array(glob.glob(join(os.environ['MD04_DIR'],version, "subhalos", "out_*_subhalos_inDistinct2.fits"))), n.array(glob.glob(join(os.environ['MD10_DIR'],version, "subhalos", "out_*_subhalos_inDistinct2.fits"))), n.array(glob.glob(join(os.environ['MD25_DIR'],version, "subhalos", "out_*_subhalos_inDistinct2.fits"))), n.array(glob.glob(join(os.environ['MD25NW_DIR'],version, "subhalos", "out_*_subhalos_inDistinct2.fits"))), n.array(glob.glob(join(os.environ['MD40_DIR'],version, "subhalos", "out_*_subhalos_inDistinct2.fits"))), n.array(glob.glob(join(os.environ['MD40NW_DIR'],version, "subhalos", "out_*_subhalos_inDistinct2.fits"))) ))
flist_3 = n.hstack(( n.array(glob.glob(join(os.environ['MD04_DIR'],version, "subhalos", "out_*_subhalos_inDistinct3.fits"))), n.array(glob.glob(join(os.environ['MD10_DIR'],version, "subhalos", "out_*_subhalos_inDistinct3.fits"))), n.array(glob.glob(join(os.environ['MD25_DIR'],version, "subhalos", "out_*_subhalos_inDistinct3.fits"))), n.array(glob.glob(join(os.environ['MD25NW_DIR'],version, "subhalos", "out_*_subhalos_inDistinct3.fits"))), n.array(glob.glob(join(os.environ['MD40_DIR'],version, "subhalos", "out_*_subhalos_inDistinct3.fits"))), n.array(glob.glob(join(os.environ['MD40NW_DIR'],version, "subhalos", "out_*_subhalos_inDistinct3.fits"))) ))

flist_1.sort()
flist_2.sort()
flist_3.sort()

def get_ids(hdB, mmin=14.5, mmax=15.5):
	msel = (hdB['mvir_cen']>mmin) & (hdB['mvir_cen']<mmax)
	return set(hdB['id_cen'][msel])
	

exponent = 4.
fsat_unev = lambda xi, a, b, N0 :  N0 * xi**a * n.e**(-b*xi**3.)
fsat = lambda xi, a, b, N0, exponent :  N0 * xi**a * n.e**(-b*xi**exponent)
logfsat= lambda logxi, a, b, logN0, exponent : n.log10( 10**logN0 * (10**logxi)**a * n.e**(-b*(10**logxi)**exponent))

def get_hist_MR(hdB, Msat = 'mvir_sat', mmin=14.5, mmax=15.5, Lbox=400.,dlogBins = 0.05, MP = 9, stat=False):
	"""return  dNsat / volume / dln(Msub/Mdistinct)
	"""
	#print hdB['mvir_cen']
	#print mmin, mmax
	msel = (hdB['mvir_cen']>mmin) & (hdB['mvir_cen']<mmax) & (hdB[Msat]>MP)
	massR = - hdB['mvir_cen'][msel] + hdB[Msat][msel]
	bins = n.arange(-6, 0.06, dlogBins)
	xb = (bins[1:]+bins[:-1])/2.
	NcenWS04 = n.histogram(massR, bins, weights=n.ones_like(massR)/Lbox**3./(dlogBins*n.log(10))*(10**(mmin/2.+mmax/2.)/rhom))[0]
	NNN,bins0 = n.histogram(massR, bins)
	#bins0 = n.histogram(massR, bins)[1]
	ok = (xb>0.3+MP-mmin)
	if stat :
		print "MD",Lbox,",Nhalo in distinct with", mmin, "<m<",mmax, "=", len(hdB['mvir_cen'][msel])
		print "bins",bins0[NNN>10]+(mmin+mmax)/2.
		print "Nsub",NNN[NNN>10] 
	return xb, NcenWS04, NNN, ok

def get_total(hdB_1, hdB_2, hdB_3, Lbox, mmin=14.5, mmax=15.5, MP=9):
	"""return  dNsat / volume / d(Msub/Mdistinct)
	print '------------------------------------------------------------------'
	print '------------------------------------------------------------------'
	"""
	#print '----------------- mvir_sat'
	xb, ratio_1, NN_1,ok_1 = get_hist_MR(hdB_1, 'mvir_sat', Lbox=Lbox, mmin=mmin, mmax=mmax, MP=MP, stat=False)
	#print '----------------- mvir_sat_sat'
	xb, ratio_2, NN_2,ok_1 = get_hist_MR(hdB_2, 'mvir_sat_n_sat_n_1', Lbox= Lbox, mmin=mmin, mmax=mmax,MP=MP, stat=False)
	#print '----------------- mvir_sat_sat_sat'
	xb, ratio_3, NN_3,ok_1 = get_hist_MR(hdB_3, 'mvir_sat_n_sat_n_1_sat_n_2', Lbox= Lbox, mmin=mmin, mmax=mmax,MP=MP, stat=False)
	
	err = (NN_1+NN_2+NN_3)**(-0.5)
	return xb, (ratio_1+ratio_2+ratio_3)*10**-xb, err, ok_1
	
def get_SHMFR(index, mmin, mmax):
	p.figure(0, (5,5))
	p.axes([0.17, 0.17, 0.75, 0.75])
	print flist_1[index]
	print flist_2[index]
	print flist_3[index]
	hdB_1 = fits.open(flist_1[index])[1].data
	hdB_2 = fits.open(flist_2[index])[1].data
	hdB_3 = fits.open(flist_3[index])[1].data
	
	boxZN = int(os.path.basename(flist_1[index]).split('_')[1])
	if flist_1[index].find('MD_0.4Gpc')>0:
		boxName='MD_0.4Gpc'
		nSN, aSN = n.loadtxt(join(os.environ['MD04_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(9.63 * 10**7)
		boxLength = 400./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		boxLengthComoving = 400.
		
	elif flist_1[index].find('MD_1Gpc')>0 :
		boxName='MD_1Gpc'
		nSN, aSN = n.loadtxt(join(os.environ['MD10_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(1.51 * 10**9)
		boxLength = 1000./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		boxLengthComoving = 1000.

	elif flist_1[index].find('MD_2.5GpcNW')>0 :
		boxName='MD_2.5GpcNW'
		nSN, aSN = n.loadtxt(join(os.environ['MD25NW_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'aSN'), 'formats': ('i4', 'f4')})
		conversion = dict(n.transpose([ nSN, 1/aSN-1 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(2.359 * 10**10)
		boxLength = 2500./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		boxLengthComoving = 2500.

	elif flist_1[index].find('MD_4GpcNW')>0 :
		boxName='MD_4GpcNW'
		nSN, redshift40, aSN = n.loadtxt(join(os.environ['MD40NW_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'redshift', 'aSN'), 'formats': ('i4', 'f4', 'f4')})
		conversion = dict(n.transpose([ nSN, redshift40 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(9.6 * 10**10)
		boxLengthComoving = 4000.
		
	elif flist_1[index].find('MD_2.5Gpc')>0 :
		boxName='MD_2.5Gpc'
		nSN, aSN, redshift25 = n.loadtxt(join(os.environ['MD25_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'aSN', 'redshift'), 'formats': ('i4', 'f4', 'f4')})
		conversion = dict(n.transpose([ nSN, redshift25 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(2.359 * 10**10)
		boxLengthComoving = 2500.
		
	elif flist_1[index].find('MD_4Gpc')>0 :
		boxName='MD_4Gpc'
		nSN, redshift40, aSN = n.loadtxt(join(os.environ['MD40_DIR'],"redshift-list.txt"), unpack=True, dtype={'names': ('nSN', 'redshift', 'aSN'), 'formats': ('i4', 'f4', 'f4')})
		conversion = dict(n.transpose([ nSN, redshift40 ]))
		boxRedshift =  conversion[boxZN] 
		logmp = n.log10(9.6 * 10**10 )
		boxLength = 4000./cosmoMD.h/cosmoMD.efunc(boxRedshift)
		boxLengthComoving = 4000.
		
	xb, y, err, ok = get_total(hdB_1, hdB_2, hdB_3, boxLengthComoving, mmin, mmax, logmp)
	ok2 = (ok)&(n.isnan(y)==False)&(y>0)&(err>0)&(err != n.inf)
	x_data = xb[ok2]
	y_data = y[ok2]
	y_data_err = err[ok2]
	
	arrr = n.ones_like(x_data)
	return n.transpose([x_data, y_data, y_data_err, boxRedshift*arrr, logmp*arrr, boxLengthComoving*arrr, mmin*arrr, mmax*arrr])
	#n.savetxt(join(os.environ['MVIR_DIR'], 'shmfr_'+str(mmin)+"_M_"+str(mmax)+".txt"), n.transpose([x_data[pouet], n.log10(y_data[pouet]), 0.05+y_data_err[pouet]]))


mmin=12.5
mmax=13.5
ttt = get_SHMFR(0, mmin, mmax) 
for ii in range(1, len(flist_1), 1):
	new = get_SHMFR(ii, 12.5, 13.5)
	print new.shape	
	if len(new)>0:
		ttt=n.vstack((ttt,new))

		
#n.ones_like(massR)/Lbox**3./(dlogBins*n.log(10))*(10**(mmin/2.+mmax/2.)/rhom))[0]
meanRHO = rhom
dlb = 0.05*n.log(10)
#1./(pow(boxLengthComoving, 3.) * pow(10, (mmin/2.+mmax/2.)))

n.savetxt(join(os.environ['MVIR_DIR'], 'shmfr_zevol_'+str(mmin)+"_M_"+str(mmax)+".txt"), ttt, header='x_data y_data y_data_err boxRedshift logmp boxLengthComoving mmin mmax')
sys.exit()

outs = []
mms = n.hstack(( n.arange(12.5, 14.6, 0.5), 15.5 ))
for mmin, mmax in zip(mms[:-1], mms[1:]):
	print mmin, mmax
	outs.append( get_SHMFR(mmin, mmax) )

for out in outs:
	print n.round(out[0][0],4), n.round(out[1].diagonal()[0]**0.5,4)

for out in outs:
	print n.round(out[0][1],4), n.round(out[1].diagonal()[1]**0.5,4)

for out in outs:
	print n.round(out[0][2],4), n.round(out[1].diagonal()[2]**0.5,4)

for out in outs:
	print n.round(out[0][3],4), n.round(out[1].diagonal()[3]**0.5,4)

import glob

datalist=n.array(glob.glob(join(os.environ['MVIR_DIR'], "shmfr_*_M_*.txt")))
x_fit=[]
y_fit=[]
yerr_fit=[]
for file in datalist:
	xx, yy, ye = n.loadtxt(file, unpack = True)
	x_fit.append(xx)
	y_fit.append(yy)
	yerr_fit.append(ye)
	
out = curve_fit(logfsat, n.hstack((x_fit)), n.hstack((y_fit)), sigma = n.hstack((yerr_fit)), p0 = p_init, maxfev = 500000000)

print out[0], out[1].diagonal()**0.5