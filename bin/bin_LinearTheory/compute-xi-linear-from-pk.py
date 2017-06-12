import astropy.units as uu
import astropy.cosmology as co
aa = co.Planck13
from scipy.interpolate import interp1d
import numpy as n
import matplotlib
#matplotlib.use('pdf')
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
from hankel import SphericalHankelTransform
	
dir = ".." #join("D:","data","MultiDark")
zList_file =  join(dir, "z-list-all-boxes.txt") 
z0 = n.loadtxt(zList_file,unpack=True)
interp1d(n.arange(len(z0)), z0 )
pklist = n.array(glob.glob(join(dir, "Pk_DM_CLASS", "MD_z*_pk.dat")))

for j, pkfile in enumerate(pklist):
	print pkfile
	snapNum = pkfile.split('_')[-2][1:]
	zSnap = z0[ int(snapNum) - 1 ]
	ks_lin_0, pks_lin_0 = n.loadtxt(pkfile, unpack=True)
	sel = (ks_lin_0<1e4)
	ks_lin, pks_lin = ks_lin_0[sel], pks_lin_0[sel]
	pk = interp1d(n.hstack((1e-20,ks_lin.min()/2., ks_lin, ks_lin.max()*2.,1e20)), n.hstack((0.,0.,pks_lin,0.,0.)))
	h = SphericalHankelTransform(nu=0,N=10000,h=0.00001)  #Create the HankelTransform instance
	#h.transform(pk)
	fr = lambda x, r : (x/r)**2. * pk(x/r)
	Rs= n.arange(0.1,200,0.1)
	xiR = n.empty_like(Rs)
	for i, R in enumerate(Rs):
		#print R, time.time()
		f = lambda x : fr(x,R)
		xiR[i] = h.transform(f)[0] / (2*n.pi**2 * R)

	n.savetxt(join(dir,"Pk_DM_CLASS","MD_z"+snapNum+"_xi.dat"),n.transpose([Rs,xiR]), header =" r xi")

