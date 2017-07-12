import astropy.io.fits as fits
import matplotlib.pyplot as p
import numpy as n
from os.path import join
import os
import sys 

summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data     

def get_files(el):
    flist_1 = join(os.environ['MD10'], "substructure", "out_"+el['snap_name']+"_subH_inDistinct_d1.fits")
    flist_2 = join(os.environ['MD10'], "substructure", "out_"+el['snap_name']+"_subH_inDistinct_d2.fits")
    flist_3 = join(os.environ['MD10'], "substructure", "out_"+el['snap_name']+"_subH_inDistinct_d3.fits")
    return flist_1, flist_2, flist_3

from hmf import MassFunction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
sigma_val=0.8229
delta_c = 1.686
volume = (1000/0.6777)**3.

from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p


index = 0
f1, f2, f3 = get_files(summ[index])
boxRedshift = summ['redshift'][index]

omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)
hf = MassFunction(cosmo_model=cosmo, sigma_8=sigma_val, z=boxRedshift, delta_h=DeltaVir_bn98(boxRedshift), delta_wrt='mean', Mmin=7, Mmax=16.5)
hf1 = MassFunction(cosmo_model=cosmo, sigma_8=sigma_val, z=1., delta_h=DeltaVir_bn98(boxRedshift), delta_wrt='mean', Mmin=7, Mmax=16.5)

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
# in units (Msun/h) / (Mpc/h)**3

ftC16 = f_BH(hf.sigma[100:-100], 0.279, 0.908, 0.671, 1.737)
rhom_units = cosmo.Om(boxRedshift)*cosmo.critical_density(boxRedshift).to(u.solMass/(u.Mpc)**3.)#/(cosmo.h)**2.
rhom = rhom_units.value # hf.mean_density#/(hz)**2. 



MF_MD = interp1d(mass, ftC16*rhom*abs(dlnsigmadlnm)/mass)

NpartMin = 50.
p_init = (-1.85, 7., -2.3, 4.)

mp04 = n.log10(NpartMin*9.63 * 10**7)
mp10 = n.log10(NpartMin*1.51 * 10**9)
mp25 = n.log10(NpartMin*2.359 * 10**10)
mp40 = n.log10(NpartMin*9.6 * 10**10. )

# fitting function
exponent = 4.
fsat_unev = lambda xi, a, b, N0 :  N0 * xi**a * n.e**(-b*xi**3.)
fsat = lambda xi, a, b, N0, exponent :  N0 * xi**a * n.e**(-b*xi**exponent)
logfsat= lambda logxi, a, b, logN0, exponent : n.log10( 10**logN0 * (10**logxi)**a * n.e**(-b*(10**logxi)**exponent))


def get_ids(hdB, mmin=14.5, mmax=15.5):
	msel = (hdB['mvir_cen']>mmin) & (hdB['mvir_cen']<mmax)
	return set(hdB['id_cen'][msel])
	

def get_hist_MR(hdB, Msat = 'mvir_sat', mmin=13.5, mmax=15.5, N_rvir=0.1, Lbox=1000.,dlogBins = 0.05, MP = 9, stat=False, rhom=1.):
	"""return  dNsat / volume / dln(Msub/Mdistinct)
	"""
	print hdB['mvir_cen']
	#print mmin, mmax
	dx2 = (hdB['x_cen']-hdB['x'+Msat[4:]])**2.
	dy2 = (hdB['y_cen']-hdB['y'+Msat[4:]])**2.
	dz2 = (hdB['z_cen']-hdB['z'+Msat[4:]])**2.
	rvir_distance = (hdB['rvir_cen'])*N_rvir/1000.#+hdB['rvir'+Msat[4:]])*N_rvir/1000.
	#print rvir_distance, (dx2+dy2+dz2)**0.5
	close_cs = ((dx2+dy2+dz2)**0.5 < rvir_distance)
	msel = (close_cs)&(hdB['mvir_cen']>mmin) & (hdB['mvir_cen']<mmax) & (hdB[Msat]>MP)
	massR = - hdB['mvir_cen'][msel] + hdB[Msat][msel]
	bins = n.arange(-6, 0.06, dlogBins)
	xb = (bins[1:]+bins[:-1])/2.
	print "rho=",rhom, ', L box=',Lbox, ', dlogBins=',dlogBins
	NcenWS04 = n.histogram(massR, bins, weights=n.ones_like(massR)/Lbox**3./(dlogBins*n.log(10))*(10**(mmin/2.+mmax/2.)/rhom))[0]
	NNN,bins0 = n.histogram(massR, bins)
	#bins0 = n.histogram(massR, bins)[1]
	ok = (xb>0.3+MP-mmin)
	if stat :
		print "in MD",Lbox,", N(halo in distinct with", mmin, "<m<",mmax, ")=", len(hdB['mvir_cen'][msel])
		print "bins",bins0[NNN>10]+(mmin+mmax)/2.
		print "Nsub",NNN[NNN>10] 
	return xb, NcenWS04, NNN, ok

def get_total(hdB_1, hdB_2, hdB_3, Lbox, mmin=14.5, mmax=15.5, MP=9, N_rvir=0.1, rhom=1.):
	"""return  dNsat / volume / d(Msub/Mdistinct)
	print '------------------------------------------------------------------'
	print '------------------------------------------------------------------'
	"""
	#print '----------------- mvir_sat'
	xb, ratio_1, NN_1,ok_1 = get_hist_MR(hdB_1, Msat='mvir_sat', mmin=mmin, mmax=mmax, N_rvir=N_rvir, Lbox=Lbox, MP=MP, stat=False, rhom=rhom)
	#print '----------------- mvir_sat_sat'
	xb, ratio_2, NN_2,ok_1 = get_hist_MR(hdB_2, Msat='mvir_sat_n_sat_n_1', mmin=mmin, mmax=mmax, N_rvir=N_rvir, Lbox=Lbox, MP=MP, stat=False, rhom=rhom)
	#print '----------------- mvir_sat_sat_sat'
	xb, ratio_3, NN_3,ok_1 = get_hist_MR(hdB_3, Msat='mvir_sat_n_sat_n_1_sat_n_2',mmin=mmin, mmax=mmax, N_rvir=N_rvir, Lbox=Lbox, MP=MP, stat=False, rhom=rhom)
	
	err = (NN_1+NN_2+NN_3)**(-0.5)
	return xb, (ratio_1+ratio_2+ratio_3)*10**-xb, err, ok_1



	
def get_SHMFR(path_f1, path_f2, path_f3, redshift, mmin, mmax, N_rvir=1., logmp = n.log10(1.51 * 10**9), boxLengthComoving = 1000., rhom=1.):
	hdB_1 = fits.open(path_f1)[1].data
	hdB_2 = fits.open(path_f2)[1].data
	hdB_3 = fits.open(path_f3)[1].data

	boxLength = 1000. # boxLengthComoving/cosmo.h/cosmo.efunc(redshift)

	xb, y, err, ok = get_total(hdB_1, hdB_2, hdB_3, Lbox=boxLengthComoving, mmin=mmin, mmax=mmax, MP=logmp, N_rvir=N_rvir, rhom=rhom)
	
	ok2 = (ok)&(n.isnan(y)==False)&(y>0)&(err>0)&(err != n.inf)
	x_data = xb[ok2]
	y_data = y[ok2]
	y_data_err = err[ok2]
	print "z=", redshift
	arrr = n.ones_like(x_data)
	return n.transpose([x_data, y_data, y_data_err, boxRedshift*arrr, logmp*arrr, boxLengthComoving*arrr, mmin*arrr, mmax*arrr])
	#n.savetxt(join(os.environ['MVIR_DIR'], 'shmfr_'+str(mmin)+"_M_"+str(mmax)+".txt"), n.transpose([x_data[pouet], n.log10(y_data[pouet]), 0.05+y_data_err[pouet]]))


index = 0
f1, f2, f3 = get_files(summ[index])
boxRedshift = summ['redshift'][index]
rhom_units = cosmo.Om(boxRedshift)*cosmo.critical_density(boxRedshift).to(u.solMass/(u.Mpc)**3.)#/(cosmo.h)**2.
rhom = rhom_units.value # hf.mean_density#/(hz)**2. 
data = get_SHMFR(f1, f2, f3, boxRedshift, 13.5, 14.5, N_rvir=1., rhom=rhom)

for el in summ[1:]:
	f1, f2, f3 = get_files(el)
	boxRedshift = el['redshift']
	# fixed to z=0 rho_m
	rhom_units = cosmo.Om(0.)*cosmo.critical_density(0.).to(u.solMass/(u.Mpc)**3.)#/(cosmo.h)**2.
	rhom = rhom_units.value # hf.mean_density#/(hz)**2. 
	# choose redshift dependent halo mass bin
	define Mstar from delta_c = sigma => halo mass bin
	0.9 Mstar < M < 1.1 Mstar
	0.09 Mstar < M < 0.11 Mstar
	
	9.9 Mstar < M < 10.1 Mstar
	
	# evolving rho_m
	#rhom_units = cosmo.Om(boxRedshift)*cosmo.critical_density(boxRedshift).to(u.solMass/(u.Mpc)**3.)#/(cosmo.h)**2.
	#rhom = rhom_units.value # hf.mean_density#/(hz)**2. 
	
	new = get_SHMFR(f1, f2, f3, boxRedshift, 13.5, 14.5, N_rvir=1., rhom=rhom)
        #print new.shape 
        if len(new)>0:
                data=n.vstack((data,new))

n.savetxt(os.path.join(os.environ['MD10'], 'substructure', 'substructure_zevol_135_Mcen_145_MD10.txt'), data, header='x_data y_data y_data_err boxRedshift logmp boxLengthComoving mmin mmax')
# log10(y_data) = log10( N/(volume . 0.05 ln(10) . M_distinct . rho_m) . M_sub/M_distinct )
# x_data = log10(M_sub/M_distinct)
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
