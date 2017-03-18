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

hd04_1 = fits.open(join(os.environ['MD04_DIR'],version, "subhalos", "out_88_subhalos_inDistinct.fits"))[1].data
hd04_2 = fits.open(join(os.environ['MD04_DIR'],version, "subhalos", "out_88_subhalos_inDistinct2.fits"))[1].data
hd04_3 = fits.open(join(os.environ['MD04_DIR'],version, "subhalos", "out_88_subhalos_inDistinct3.fits"))[1].data
mp04 = n.log10(NpartMin*9.63 * 10**7)

hd10_1 = fits.open(join(os.environ['MD10_DIR'],version, "subhalos", "out_128_subhalos_inDistinct.fits"))[1].data
hd10_2 = fits.open(join(os.environ['MD10_DIR'],version, "subhalos", "out_128_subhalos_inDistinct2.fits"))[1].data
hd10_3 = fits.open(join(os.environ['MD10_DIR'],version, "subhalos", "out_128_subhalos_inDistinct3.fits"))[1].data
mp10 = n.log10(NpartMin*1.51 * 10**9)

hd25_1 = fits.open(join(os.environ['MD25_DIR'],version, "subhalos", "out_80_subhalos_inDistinct.fits"))[1].data
hd25_2 = fits.open(join(os.environ['MD25_DIR'],version, "subhalos", "out_80_subhalos_inDistinct2.fits"))[1].data
hd25_3 = fits.open(join(os.environ['MD25_DIR'],version, "subhalos", "out_80_subhalos_inDistinct3.fits"))[1].data
mp25 = n.log10(NpartMin*2.359 * 10**10)

hd25nw_1 = fits.open(join(os.environ['MD25NW_DIR'],version, "subhalos", "out_80_subhalos_inDistinct.fits"))[1].data
hd25nw_2 = fits.open(join(os.environ['MD25NW_DIR'],version, "subhalos", "out_80_subhalos_inDistinct2.fits"))[1].data
hd25nw_3 = fits.open(join(os.environ['MD25NW_DIR'],version, "subhalos", "out_80_subhalos_inDistinct3.fits"))[1].data
mp25nw = mp25

hd40_1 = fits.open(join(os.environ['MD40_DIR'],version, "subhalos", "out_128_subhalos_inDistinct.fits"))[1].data
hd40_2 = fits.open(join(os.environ['MD40_DIR'],version, "subhalos", "out_128_subhalos_inDistinct2.fits"))[1].data
hd40_3 = fits.open(join(os.environ['MD40_DIR'],version, "subhalos", "out_128_subhalos_inDistinct3.fits"))[1].data
mp40 = n.log10(NpartMin*9.6 * 10**10. )

hd40nw_1 = fits.open(join(os.environ['MD40NW_DIR'],version, "subhalos", "out_16_subhalos_inDistinct.fits"))[1].data
hd40nw_2 = fits.open(join(os.environ['MD40NW_DIR'],version, "subhalos", "out_16_subhalos_inDistinct2.fits"))[1].data
hd40nw_3 = fits.open(join(os.environ['MD40NW_DIR'],version, "subhalos", "out_16_subhalos_inDistinct3.fits"))[1].data
mp40nw = mp40


def get_ids(hd04_1, mmin=14.5, mmax=15.5):
	msel = (hd04_1['mvir_cen']>mmin) & (hd04_1['mvir_cen']<mmax)
	return set(hd04_1['id_cen'][msel])
	
#id_1=get_ids(hd04_1)
#id_2=get_ids(hd04_2)
#id_3=get_ids(hd04_3)

#hd04_1['GroupSize'][msel]
#hd04_1['GroupID'][msel]

allidsat = set(hd04_1['id_sat'])

fsat_unev = lambda xi, a, b, N0 :  N0 * xi**a * n.e**(-b*xi**3.)
fsat = lambda xi, a, b, N0 :  N0 * xi**a * n.e**(-b*xi**4.)
logfsat= lambda logxi, a, b, logN0 : n.log10( 10**logN0 * (10**logxi)**a * n.e**(-b*(10**logxi)**4.))

def get_hist_MR(hd04_1, Msat = 'mvir_sat', mmin=14.5, mmax=15.5, Lbox=400.,dlogBins = 1., MP = 9, stat=False):
	"""return  dNsat / volume / dln(Msub/Mdistinct)
	"""
	msel = (hd04_1['mvir_cen']>mmin) & (hd04_1['mvir_cen']<mmax) & (hd04_1[Msat]>MP)
	massR = - hd04_1['mvir_cen'][msel] + hd04_1[Msat][msel]
	bins = n.arange(-6, 0.06, dlogBins)
	xb = (bins[1:]+bins[:-1])/2.
	NcenWS04 = n.histogram(massR, bins, weights=n.ones_like(massR)/Lbox**3./(dlogBins*n.log(10))*(10**(mmin/2.+mmax/2.)/rhom))[0]
	NNN,bins0 = n.histogram(massR, bins)
	#bins0 = n.histogram(massR, bins)[1]
	ok = (xb>0.3+MP-mmin)
	if stat :
		print "MD",Lbox,",Nhalo in distinct with", mmin, "<m<",mmax, "=", len(hd04_1['mvir_cen'][msel])
		print "bins",bins0[NNN>10]+(mmin+mmax)/2.
		print "Nsub",NNN[NNN>10] 
	return xb, NcenWS04, NNN, ok

def get_total(hd04_1, hd04_2, hd04_3, Lbox, mmin=14.5, mmax=15.5, MP=9):
	"""return  dNsat / volume / d(Msub/Mdistinct)
	print '------------------------------------------------------------------'
	print '------------------------------------------------------------------'
	"""
	print '----------------- mvir_sat'
	xb, ratio_1, NN_1,ok_1 = get_hist_MR(hd04_1, 'mvir_sat', Lbox=Lbox, mmin=mmin, mmax=mmax, MP=MP, stat=True)
	print '----------------- mvir_sat_sat'
	xb, ratio_2, NN_2,ok_1 = get_hist_MR(hd04_2, 'mvir_sat_n_sat_n_1', Lbox= Lbox, mmin=mmin, mmax=mmax,MP=MP, stat=True)
	print '----------------- mvir_sat_sat_sat'
	xb, ratio_3, NN_3,ok_1 = get_hist_MR(hd04_3, 'mvir_sat_n_sat_n_1_sat_n_2', Lbox= Lbox, mmin=mmin, mmax=mmax,MP=MP, stat=True)
	
	err = (NN_1+NN_2+NN_3)**(-0.5)
	return xb, (ratio_1+ratio_2+ratio_3)*10**-xb, err, ok_1
	
def plot_SHMFR(mmin, mmax):
	p.figure(0, (5,5))
	p.axes([0.17, 0.17, 0.75, 0.75])
	print '------------------------------------------------------------------'
	print 'MD04'
	print '------------------------------------------------------------------'
	xb, y, err, ok = get_total(hd04_1, hd04_2, hd04_3, 400., mmin, mmax, mp04)
	#print ok
	x_data = xb[ok]
	y_data = y[ok]
	y_data_err = err[ok]
	if len(xb[ok])>2:
		#print len(xb[ok])
		p.errorbar(xb[ok], n.log10(y[ok])+xb[ok], yerr= err[ok], label='M04')
	
	print '------------------------------------------------------------------'
	print 'MD10'
	print '------------------------------------------------------------------'
	xb, y, err, ok = get_total(hd10_1, hd10_2, hd10_3, 1000., mmin, mmax, mp10)
	#print ok
	if len(xb[ok])>2:
		p.errorbar(xb[ok], n.log10(y[ok])+xb[ok], yerr= err[ok], label='M10')
		x_data = n.hstack((x_data, xb[ok]))
		y_data = n.hstack((y_data, y[ok]))
		y_data_err = n.hstack((y_data_err, err[ok]))

		
	print '------------------------------------------------------------------'
	print 'MD25'
	print '------------------------------------------------------------------'
	xb, y, err, ok = get_total(hd25_1, hd25_2, hd25_3, 2500., mmin, mmax, mp25)
	#print ok
	if len(xb[ok])>2:
		p.errorbar(xb[ok], n.log10(y[ok])+xb[ok], yerr= err[ok], label='M25')
		x_data = n.hstack((x_data, xb[ok]))
		y_data = n.hstack((y_data, y[ok]))
		y_data_err = n.hstack((y_data_err, err[ok]))

	print '------------------------------------------------------------------'
	print 'MD25n'
	print '------------------------------------------------------------------'
	xb, y, err, ok = get_total(hd25nw_1, hd25nw_2, hd25nw_3, 2500., mmin, mmax, mp25nw)
	#print ok
	if len(xb[ok])>2:
		p.errorbar(xb[ok], n.log10(y[ok])+xb[ok], yerr= err[ok], label='M25n')
		x_data = n.hstack((x_data, xb[ok]))
		y_data = n.hstack((y_data, y[ok]))
		y_data_err = n.hstack((y_data_err, err[ok]))

	print '------------------------------------------------------------------'
	print 'MD40'
	print '------------------------------------------------------------------'
	xb, y, err, ok = get_total(hd40_1, hd40_2, hd40_3, 4000., mmin, mmax, mp40)
	#print ok
	if len(xb[ok])>2:
		p.errorbar(xb[ok], n.log10(y[ok])+xb[ok], yerr= err[ok], label='M40')
		x_data = n.hstack((x_data, xb[ok]))
		y_data = n.hstack((y_data, y[ok]))
		y_data_err = n.hstack((y_data_err, err[ok]))

	print '------------------------------------------------------------------'
	print 'MD40n'
	print '------------------------------------------------------------------'
	xb, y, err, ok = get_total(hd40nw_1, hd40nw_2, hd40nw_3, 4000., mmin, mmax, mp40nw)
	#print ok
	if len(xb[ok])>2:
		p.errorbar(xb[ok], n.log10(y[ok])+xb[ok], yerr= err[ok], label='M40n')
		x_data = n.hstack((x_data, xb[ok]))
		y_data = n.hstack((y_data, y[ok]))
		y_data_err = n.hstack((y_data_err, err[ok]))

	
	pouet = (y_data>0)
	out = curve_fit(logfsat, x_data[pouet], n.log10(y_data[pouet]), sigma = 0.05+y_data_err[pouet], p0 = (-1.85, 7., -2.3), maxfev = 500000000) 
	print out[0], out[1].diagonal()**0.5
	xx = n.arange(-6,0, 0.01)
	#p.plot(xx, n.log10(fsat_unev(10**xx, -1.8, 6.283, 0.21)/(10**(mmin/2.+mmax/2.)/rhom))+xx, label='unevolved', ls='solid', color='k')
	p.plot(xx, logfsat(xx, out[0][0], out[0][1], out[0][2])+xx, label='fit', ls='dashed', color='k')
	p.ylabel(r'$\log_{10}\left[ \frac{M_d M_s}{\rho_m} \frac{dn}{dM_s} \right] $') 
	p.xlabel(r'$\log_{10}(M_{s}/M_{d})$')
	p.title(r"$"+str(mmin)+"<M_{s}<"+str(mmax)+"$")
	p.legend(loc=0, frameon=False)
	#p.yscale('log')
	p.ylim((-5, 1))
	p.xlim(( -4, 0 )) 
	p.grid()
	#p.savefig(join(os.environ['MVIR_DIR'], 'shmfr_'+str(mmin)+"_M_"+str(mmax)+".png"))
	p.clf()
	return out

outs = []
mms = n.hstack(( n.arange(12.5, 14.6, 0.5), 15.5 ))
for mmin, mmax in zip(mms[:-1], mms[1:]):
	print mmin, mmax
	outs.append( plot_SHMFR(mmin, mmax) )

for out in outs:
	print n.round(out[0][0],4), n.round(out[1].diagonal()[0]**0.5,4)

for out in outs:
	print n.round(out[0][1],4), n.round(out[1].diagonal()[1]**0.5,4)

for out in outs:
	print n.round(out[0][2],4), n.round(out[1].diagonal()[2]**0.5,4)
