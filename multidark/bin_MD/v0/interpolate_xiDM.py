import numpy as n
from scipy.interpolate import interp1d,interp2d
from scipy.integrate import quad
from scipy.optimize import curve_fit
from scipy.special import gamma
import glob
import pylab as p
import scipy.fftpack as f
"""
8=0.5431
A0=750.*(s8/0.8)**3.75
R=26.*(s8/0.8)**0.15

R01h=1.87*(s8/0.8)**-0.47
R1h=3.87*(s8/0.8)**0.29
R1=3.33*(s8/0.8)**0.88
R2h=1.69*(s8/0.8)**0.43

xi_seljak=lambda r : -A0*n.e**(-r/R)/(4*n.pi*r*R**2.)*( 1- (R/R1h)**2*n.e**(-(R+R1h)*r/(R*R1h)))
xiBB_seljak=lambda r : -A0*n.e**(-r/R)/(4*n.pi*r*R**2.)
"""

def interpolate_xi_fft(pk,sampling=0.001):
	ks=n.arange(pk.x.min(),pk.x.max(),sampling)
	pkL=n.array([pk(ks) for k in ks])
	nn=len(pkL)
	four=nn*sampling*f.ifft(pkL*ks)/(2.*n.pi**2.)
	rrs=n.arange(nn)*2*n.pi/(nn*sampling)
	#print four.imag,rrs
	return interp1d(rrs,four.imag/rrs,kind='linear')

def interpolate_xi_windowed(pk,rxi,kmin=0.0002,kmax=40.):
	def xi_fun(r):
		pkm_integrand=lambda k,r: pk(k) * k**3./(2.*n.pi**2.) * n.sin(k*r) / (k**2.0 * r) / k
		oout=quad(pkm_integrand,kmin,kmax,args=(r),limit=20000, epsabs=1.49e-06, epsrel=1.49e-06)
		return oout[0]/(4.*n.pi)
	xis=n.array([xi_fun(r) for r in rxi])	
	return interp1d(rxi,xis)

nSnap,zSnap,aSnap=n.loadtxt("/Volumes/data/BigMD/1Gpc_3840_Planck1/redshift-snapshot.list",unpack=True,dtype={'names': ('nSnap', 'zSnap', 'aSnap'),'formats': ('i4', 'f4', 'f4')})
aSnap=aSnap[n.argsort(zSnap)]
nSnap=nSnap[n.argsort(zSnap)]
zSnap=zSnap[n.argsort(zSnap)]

#snS=[19,20,21,22,23,24,25,27,36,87]
snS=n.array([ 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 47, 49, 51, 53, 56, 60, 66, 70, 74, 76, 78, 80, 82, 87])

for sn in snS:
	#index=n.searchsorted(zSnap,0.78)
	#zSnap[index-1]
	print sn
	pkDMf=glob.glob("/Volumes/data/BigMD/1Gpc_3840_Planck1/Pk_DM/Delta2_folded_*"+ str(sn) +".txt")[0]
	kDM_f,d2kDM_f=n.loadtxt(pkDMf,unpack=True,usecols=(0,1),dtype={'names': ('k', 'd2k'),'formats': ('f4', 'f4')})
	D2kDM_f=interp1d(kDM_f,d2kDM_f)

	pkDM=glob.glob( "/Volumes/data/BigMD/1Gpc_3840_Planck1/Pk_DM/DELTA2_all*"+ str(sn) +".txt")[0]
	kDM,d2kDM=n.loadtxt(pkDM,unpack=True,usecols=(0,1), dtype={'names': ('k', 'd2k'), 'formats': ('f4', 'f4')})
	D2kDM=interp1d(kDM,d2kDM)

	kk,pkDM=n.loadtxt(glob.glob("/Applications/camb/MDPL/*"+str(sn)+ "__matterpower.dat")[0], unpack=True, usecols=(0,1), dtype={'names': ('k', 'pk'),'formats': ('f4', 'f4')})

	DkDML1=interp1d(kk,kk**3.*pkDM/(2.*n.pi**2.))

	k1,k2=0.1,0.9
	factor=n.mean(D2kDM(10**n.arange(n.log10(k1), n.log10(k2),0.005))/ DkDML1( 10**n.arange (n.log10(k1), n.log10(k2),0.005)))

	DkDML=interp1d(kk,DkDML1.y*factor)

	kmin=0.0001
	kTRL=0.15
	kTRf=15.
	kmax=45.

	kss=10**n.arange(n.log10(kmin),n.log10(kmax),0.005)

	#final interpolation of D2K of the dark matter
	D2_DM=interp1d(kss,n.hstack((DkDML(kss[(kss<kTRL)]), D2kDM(kss[(kss>=kTRL)&(kss<kTRf)]), D2kDM_f(kss[(kss>=kTRf)]))))

	# transition : 0.3<k<1
	"""
	p.loglog(D2kDM_f.x,D2kDM_f.y,label='D2kDM_f')
	p.loglog(D2kDM.x,D2kDM.y,label='D2kDM')
	p.loglog(PkDML.x,PkDML.y,label='PkDM')
	p.loglog(D2_DM.x,D2_DM.y,label="interpolation")
	p.legend()
	p.show()
	xiDMS=glob.glob("/Volumes/data/ugri_clustering/xiDM_realSpace/*.dat")
	xDM,yDM=n.loadtxt(xiDMS[-1],unpack=True)
	"""

	pk=interp1d(D2_DM.x,D2_DM.y*(2.*n.pi**2.)/D2_DM.x**3.)
	n.savetxt("/Volumes/data/BigMD/1Gpc_3840_Planck1/xi_DM/Pk_" +str(sn)+ ".txt", n.transpose([pk.x,pk.y]) ,header=" k pk ")	
	print "pk done"
	#xi5=interpolate_xi_fft(pk,sampling=0.005)
	#n.savetxt("/Volumes/data/BigMD/1Gpc_3840_Planck1/xi_DM/2PCF_realspace_" +str(sn)+ ".txt", n.transpose([xi5.x[(xi5.x>0)&(xi5.x<300)],xi5.y[0][(xi5.x>0)&(xi5.x<300)]]) ,header=" r_Mpc_per_h xi_r ")

import sys
sys.exit()


xi1=interpolate_xi_fft(pk,sampling=0.2)
xi2=interpolate_xi_fft(pk,sampling=0.1)
xi3=interpolate_xi_fft(pk,sampling=0.05)
xi4=interpolate_xi_fft(pk,sampling=0.01)

p.loglog(xi1.x,xi1.y[0],label='fft 0.2')
p.loglog(xi2.x,xi2.y[0],label='fft 0.1')
p.loglog(xi3.x,xi3.y[0],label='fft 0.05')
p.loglog(xi4.x,xi4.y[0],label='fft 0.01')
p.loglog(xi5.x,xi5.y[0],label='fft 0.005')
p.loglog(xDM,yDM,label='MD')

p.legend()
p.show()


# ks=10**n.arange(n.log10(pk.x[0]),n.log10(pk.x[-1]),0.1)
# pkds=interp1d(ks,pk(ks))


pkm_integrand=lambda lgk,r: n.log(10)* n.e**(2*lgk/n.log(10)) *pk(n.e**(lgk/n.log(10))) * n.sin(n.e**(lgk/n.log(10)) *r)
rrs=10**n.arange(-2,1,0.1)
xir=n.array([quad(pkm_integrand,-3.5,1.5,args=(rr,),limit=20000, epsabs=1.49e-06, epsrel=1.49e-06)[0]/(rr*2*n.pi**2) for rr in rrs ])

dkm_integrand=lambda lgk,r: n.log(10)* n.e**(-lgk/n.log(10)) *D2_DM(n.e**(lgk/n.log(10))) * n.sin(n.e**(lgk/n.log(10)) *r)/r
xir2=n.array([quad(dkm_integrand,-3.5,1.5,args=(rr,),limit=20000, epsabs=1.49e-06, epsrel=1.49e-06)[0] for rr in rrs ])

p.loglog(rrs,xir,label='tfpk')
p.loglog(rrs,xir2,label='tfdk')

p.loglog(xDM,yDM,label='MD')

p.legend()
p.show()


fou=f.ifft(pk.y)
out=f.ifft2([pk.x,pk.y])

p.loglog(abs(out[0]),abs(out[1]))

fou.imag

xi1=interpolate_xi_fft(pkds,sampling=1.)
p.loglog(xi1.x,xi1.x*xi1.x*xi1.y[0],label='fft 0.001')
xi1=interpolate_xi_fft(pkds,sampling=0.002)
p.loglog(xi1.x,xi1.x*xi1.x*xi1.y[0],label='fft 0.002')
xi1=interpolate_xi_fft(pkds,sampling=0.004)
p.loglog(xi1.x,xi1.x*xi1.x*xi1.y[0],label='fft 0.004')
xi1=interpolate_xi_fft(pkds,sampling=0.008)
p.loglog(xi1.x,xi1.x*xi1.x*xi1.y[0],label='fft 0.008')
xi1=interpolate_xi_fft(pk,sampling=0.01)
p.loglog(xi1.x,xi1.x*xi1.x*xi1.y[0],label='fft 0.01')

rrs,xir
p.loglog(xDM,xDM*xDM*yDM,label='MD')
pk_e=interp1d(n.hstack((kk,2)),n.hstack((pkk,pk(1)/2.)))
xi1=interpolate_xi_fft(pk_e,sampling=0.001)
p.loglog(xi1.x,xi1.x*xi1.x*xi1.y[0],label='E fft 0.005')

p.legend()
p.show()




import sys
sys.exit()


xi2=interpolate_xi_windowed(pk,xDM,kmin=0.0002,kmax=1.)

pk_e=interp1d(n.hstack((kk,20)),n.hstack((pkk,0)))
xi1ext=interpolate_xi_fft(pk_e,sampling=0.1)
xi2ext=interpolate_xi_windowed(pk_e,xDM,kmin=0.0002,kmax=10.)


p.loglog(xDM,yDM,label='MD')
p.loglog(xi1.x,xi1.y[0],label='fft')
p.loglog(xi2.x,xi2.y,label='integration')
p.loglog(xDM,xi_seljak(xDM),label='Seljak')
p.legend()
p.show()
