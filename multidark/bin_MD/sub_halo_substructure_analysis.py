import astropy.io.fits as fits
import numpy as n
from os.path import join
import os
import sys 

from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
sigma_val=0.8229
delta_c = 1.686
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use('Agg')
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
p_init = (-1.85, -2.3)

env = os.environ['MD10']
sat_in_cen_d1 =join(env, "substructure", "out_0.74980_subH_inDistinct_d1.fits")
sat_in_cen_d2 =join(env, "substructure", "out_0.74980_subH_inDistinct_d2.fits")
sat_in_cen_d3 =join(env, "substructure", "out_0.74980_subH_inDistinct_d3.fits")


hd10_1 = fits.open(sat_in_cen_d1)[1].data
hd10_2 = fits.open(sat_in_cen_d2)[1].data
hd10_3 = fits.open(sat_in_cen_d3)[1].data
mp10 = n.log10(NpartMin*1.51 * 10**9)

def get_ids(hd04_1, mmin=14.5, mmax=15.5):
	msel = (hd04_1['mvir_cen']>mmin) & (hd04_1['mvir_cen']<mmax)
	return set(hd04_1['id_cen'][msel])
	
exponent = 4.
fsat_unev = lambda xi, a, b, N0 :  N0 * xi**a * n.e**(-b*xi**3.)
fsat = lambda xi, a, b, N0, exponent :  N0 * xi**a * n.e**(-b*xi**exponent)
logfsat= lambda logxi, a, logN0 : n.log10( 10**logN0 * (10**logxi)**a )#* n.e**(-b*(10**logxi)**exponent))

def get_hist_MR(hd04_1, Msat = 'mvir_sat', mmin=14.5, mmax=15.5, Lbox=400.,dlogBins = 0.05, MP = 9, stat=False):
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
	return xb, (ratio_1+ratio_2+ratio_3)*10**-xb, err, ok_1, ratio_1*10**-xb, ratio_2*10**-xb, ratio_3*10**-xb
	
def plot_SHMFR(mmin, mmax):
	p.figure(0, (5,5))
	p.axes([0.17, 0.17, 0.75, 0.75])
	print '------------------------------------------------------------------'
	print 'MD10'
	print '------------------------------------------------------------------'
	xb, y, err, ok, y_1, y_2, y_3 = get_total(hd10_1, hd10_2, hd10_3, 1000., mmin, mmax, mp10)
	print ok
	x_data = xb[ok]
	y_data = y[ok]
	y_data_1 = y_1[ok]
	y_data_2 = y_2[ok]
	y_data_3 = y_3[ok]
	y_data_err = err[ok]
	if len(xb[ok])>2:
		#print len(xb[ok])
		#p.errorbar(xb[ok], n.log10(y[ok])+xb[ok], yerr= err[ok], label='all')
		p.errorbar(xb[ok], n.log10(y_1[ok])+xb[ok], yerr= err[ok], label='d1')
		p.errorbar(xb[ok], n.log10(y_2[ok])+xb[ok], yerr= err[ok], label='d2')
		p.errorbar(xb[ok], n.log10(y_3[ok])+xb[ok], yerr= err[ok], label='d3')
	
	xx = n.arange(-6,0, 0.01)
	
	pouet = (y_data_1>0)
	if len(x_data[pouet])>10:
		out = curve_fit(logfsat, x_data[pouet], n.log10(y_data_1[pouet]), sigma = 0.05+y_data_err[pouet], p0 = p_init, maxfev = 500000000) 
		print "fit:", out[0], out[1].diagonal()**0.5
		p.plot(xx, logfsat(xx, out[0][0], out[0][1])+xx, label='fit d1 '+str(n.round(out[0][0]),2), ls='dashed', color='k')
		
	pouet = (y_data_2>0)
	if len(x_data[pouet])>10:
		out = curve_fit(logfsat, x_data[pouet], n.log10(y_data_2[pouet]), sigma = 0.05+y_data_err[pouet], p0 = p_init, maxfev = 500000000) 
		print "fit:", out[0], out[1].diagonal()**0.5
		p.plot(xx, logfsat(xx, out[0][0], out[0][1])+xx, label='fit d2 '+str(n.round(out[0][0]),2), ls='dashed', color='k')
	
	pouet = (y_data_3>0)
	if len(x_data[pouet])>10:
		out = curve_fit(logfsat, x_data[pouet], n.log10(y_data_3[pouet]), sigma = 0.05+y_data_err[pouet], p0 = p_init, maxfev = 500000000) 
		print "fit:", out[0], out[1].diagonal()**0.5
		p.plot(xx, logfsat(xx, out[0][0], out[0][1])+xx, label='fit d3 '+str(n.round(out[0][0]),2), ls='dashed', color='k')
	
	
	p.ylabel(r'$\log_{10}\left[ \frac{M_d M_s}{\rho_m} \frac{dn}{dM_s} \right] $') 
	p.xlabel(r'$\log_{10}(M_{s}/M_{d})$')
	p.title(r"$"+str(mmin)+"<M_{d}<"+str(mmax)+"$")
	p.legend(loc=0, frameon=False)
	#p.yscale('log')
	p.ylim((-7, -2))
	p.xlim(( -3, 0 )) 
	p.grid()
	p.savefig(join(os.environ['MD10'], 'substructure', 'shmfr_'+str(mmin)+"_M_"+str(mmax)+".png"))
	n.savetxt(join(os.environ['MD10'], 'substructure', 'shmfr_'+str(mmin)+"_M_"+str(mmax)+".txt"), n.transpose([x_data[pouet], n.log10(y_data[pouet]), 0.05+y_data_err[pouet]]))
	p.clf()
	return out

outs = []
mms = n.arange(13.5, 14.6, 0.5)
for mmin, mmax in zip(mms[:-1], mms[1:]):
	print mmin, mmax
	outs.append( plot_SHMFR(mmin, mmax) )

for out in outs:
	print n.round(out[0][0],4), n.round(out[1].diagonal()[0]**0.5,4)

for out in outs:
	print n.round(out[0][1],4), n.round(out[1].diagonal()[1]**0.5,4)

for out in outs:
	print n.round(out[0][2],4), n.round(out[1].diagonal()[2]**0.5,4)

for out in outs:
	print n.round(out[0][3],4), n.round(out[1].diagonal()[3]**0.5,4)

import glob

datalist=n.array(glob.glob(join(os.environ['MD10'], 'substructure', "shmfr_*_M_*.txt")))
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