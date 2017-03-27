import astropy.units as uu
from scipy.interpolate import interp1d
import numpy as n
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
from hankel import SphericalHankelTransform
import os
import lib_functions_1pt as lib

version = 'v4'
qty = "clustering"
afactor=1.
redshift = 0.
Rs, xiR  = n.loadtxt(join(os.environ['MD_DIR'],"Pk_DM_CLASS","MD_z23_xi_2017.dat"), unpack = True)
linXi = interp1d(Rs,xiR)
#print Rs, xiR		
fileXi = n.array(glob.glob( join(os.environ['MD_DIR'], "MD_*Gpc*",  version, qty,"out_*30_xiR.pkl")))
fileXi.sort()

NminParticles = 1000.
rmin = 8.
rmax=20.

print fileXi

def getBias(file):
	fileN = os.path.basename(file)
	boxN = file.split('\\')[3]
	boxZN = fileN.split('_')[1]
	f=open(file,'r')
	bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbinsL, vbinsH = cPickle.load(f)
	f.close()
	
	rr = (bin_xi3D[1:] + bin_xi3D[:-1])/2.
	# print file
	m0=float(fileN.split('_')[3])
	m1=float(fileN.split('_')[4])
	
	m_low= 10**(m0)/lib.cosmoMD.h
	m_high = 10**(m1)/lib.cosmoMD.h
	m_mean = 10**((m0 + m1)/2.)/lib.cosmoMD.h
	
	hf, boxLength, boxName, boxRedshift, logmp, boxLengthComoving, massCorrection = lib.get_basic_info(file, int(boxZN), delta_wrt='mean')
	print logmp
	if m_mean>NminParticles*10**logmp :
		m2sigma = interp1d(hf.M, hf.sigma )
		sig_low = m2sigma( m_low)
		sig_high = m2sigma( m_high)
		sig_mean = m2sigma( m_mean)
		
		xi = DR*volume/(dV * pairCount) -1.
		ok = (rr>rmin)&(rr<rmax)
		scale = (n.min(rr[ok])+n.max(rr[ok]))/2.
		bias = n.mean(xi[ok]/linXi(rr[ok]))
		biasErr = n.std(xi[ok]/linXi(rr[ok]))
		#print [sig_low, sig_high, sig_mean, scale, bias, biasErr, volume, afactor]
		print bias
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

n.savetxt(join(os.environ['MVIR_DIR'],  "halo-bias-measurement-summary.data"), tableOut, header=" sigma_low sigma_high sigma_mean scale bias biasErr volume a logmp")
