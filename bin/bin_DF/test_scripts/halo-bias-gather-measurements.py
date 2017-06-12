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
 
dir = ".." #join("D:","data","MultiDark")
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")

zList_file =  join(dir_25, "redshift_list.txt") 
n0,z0,a0 = n.loadtxt(zList_file,unpack=True)
n1 = n0.astype(int)+1
n2 = n1.astype(str)
snap2Z = {n2[i].zfill(3): z0[i] for i in range(len(n2))}
a2snap = {a0[i] : n2[i] for i in range(len(n2))}

vmins = n.array(["65.0", "72.9", "81.8", "91.8", "103.0", "115.5", "129.6", "145.5", "163.2", "183.1", "205.5", "230.6","258.7", "290.3", "325", "365.5", "410.1", "460.1", "516.3", "579.3", "650.0", "729.3", "818.3", "918.1", "1030.1", "1155.8", "1296.9", "1455.17"])
vlow = vmins.astype('float')[:-1]
vup = vmins.astype('float')[1:]
vmean = (vup*vlow)**0.5

limits_04 = [70, 2000]
limits_10 = [200, 5000]
limits_25 = [300, 5000]
limits_40 = [500, 5000]

list44 = n.array(glob.glob(join("..","MD_1Gpc","halo_bias","clustering","hlist_*_vmax_*rmax_050_xiR.pkl")))
list42 = n.array(glob.glob(join("..","MD_0.4Gpc","halo_bias","clustering","hlist_*_vmax_*rmax_015_xiR.pkl")))
list41 = n.array(glob.glob(join("..","MD_2.5Gpc","halo_bias","clustering","hlist_*_vmax_*rmax_50_xiR.pkl")))
list43 = n.array(glob.glob(join("..","MD_4Gpc","halo_bias","clustering","hlist_*_vmax_*rmax_050_xiR.pkl")))
list40=n.hstack((list41, list42, list43, list44))

def getBias(file):
		f=open(file,'r')
		bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbinsL, vbinsH = cPickle.load(f)
		f.close()
		rr = (bin_xi3D[1:] + bin_xi3D[:-1])/2.
		# print file
		aon = float(file.split('_')[-7])
		vlow= float(file.split('_')[-5])
		vhigh = float(file.split('_')[-4])
		vmean = (vlow * vhigh)**0.5
		if file.split('\\')[1] == "MD_0.4Gpc" and vlow>limits_04[0]:
			volume = 400**3.
			xi = DR*volume/(dV * pairCount) -1.
			ok = (rr>8)&(rr<20)
			vol = volume
			Rs, xiR  = n.loadtxt(join(dir,"Pk_DM_CLASS","multidark00_z"+a2snap[a0[n.searchsorted(a0,aon)]]+"_xi.dat"), unpack = True)
			linXi = interp1d(Rs,xiR)
			scale = (n.min(rr[rr>8])+n.max(rr[rr>8]))/2.
			bias = n.mean(xi[(rr>8)]/linXi(rr[(rr>8)]))
			biasErr = n.std(xi[(rr>8)]/linXi(rr[(rr>8)]))
			# print [vlow, vhigh,vmean, scale, bias, biasErr, vol, aon]
			return [vlow, vhigh,vmean, scale, bias, biasErr, vol, aon]
			
		if file.split('\\')[1] == "MD_1Gpc" and vlow>limits_10[0]:
			volume = 1000**3.
			xi = DR*volume/(dV * pairCount) -1.
			ok = (rr>8)&(rr<20)
			vol = volume
			Rs, xiR  = n.loadtxt(join(dir,"Pk_DM_CLASS","multidark00_z"+a2snap[a0[n.searchsorted(a0,aon)]]+"_xi.dat"), unpack = True)
			linXi = interp1d(Rs,xiR)
			scale = (n.min(rr[rr>8])+n.max(rr[rr>8]))/2.
			bias = n.mean(xi[(rr>8)]/linXi(rr[(rr>8)]))
			biasErr = n.std(xi[(rr>8)]/linXi(rr[(rr>8)]))
			# print [vlow, vhigh,vmean, scale, bias, biasErr, vol, aon]
			return [vlow, vhigh,vmean, scale, bias, biasErr, vol, aon]
			
		if file.split('\\')[1] == "MD_2.5Gpc" and vlow>limits_25[0]:
			n25, a25 = n.loadtxt(join("..","MD_2.5Gpc","redshift-list.txt"),unpack=True)
			aa = a25[(n25 == aon)][0]
			volume = 2500**3.
			xi = DR*volume/(dV * pairCount) -1.
			ok = (rr>8)&(rr<20)
			vol = volume
			Rs, xiR  = n.loadtxt(join(dir,"Pk_DM_CLASS","multidark00_z"+a2snap[a0[n.searchsorted(a0,aa)]]+"_xi.dat"), unpack = True)
			linXi = interp1d(Rs,xiR)
			scale = (n.min(rr[rr>8])+n.max(rr[rr>8]))/2.
			bias = n.mean(xi[(rr>8)]/linXi(rr[(rr>8)]))
			biasErr = n.std(xi[(rr>8)]/linXi(rr[(rr>8)]))
			# print [vlow, vhigh,vmean, scale, bias, biasErr, vol, a25[(n25 == aon)][0] ]
			return [vlow, vhigh,vmean, scale, bias, biasErr, vol,  aa]
			
		if file.split('\\')[1] == "MD_4Gpc" and vlow>limits_40[0]:
			n40, a40 = n.loadtxt(join("..","MD_4Gpc","redshift-list.txt"),unpack=True)
			aa = a40[(n40 == aon)][0]
			volume = 4000**3.
			xi = DR*volume/(dV * pairCount) -1.
			ok = (rr>8)&(rr<20)
			vol = volume
			Rs, xiR  = n.loadtxt(join(dir,"Pk_DM_CLASS","multidark00_z"+a2snap[a0[n.searchsorted(a0,aa)]]+"_xi.dat"), unpack = True)
			linXi = interp1d(Rs,xiR)
			scale = (n.min(rr[rr>8])+n.max(rr[rr>8]))/2.
			bias = n.mean(xi[(rr>8)]/linXi(rr[(rr>8)]))
			biasErr = n.std(xi[(rr>8)]/linXi(rr[(rr>8)]))
			# print [vlow, vhigh,vmean, scale, bias, biasErr, vol, a40[(n40 == aon)][0] ]
			return [vlow, vhigh,vmean, scale, bias, biasErr, vol, aa ]

out=[]
for el in list40 : #range(len(list40)):
	got = getBias(el)
	if got == None :
		pass
	else:	
		out.append(got)

out= n.array(out)
vlow, vhigh,vmean, scale, bias, biasErr, vol, aon = n.transpose(out)

n.savetxt(join("..", "clustering", "halo-bias-measurement-summary.data"), n.transpose([vlow, vhigh,vmean, scale, bias**0.5, 0.5*biasErr*bias**(-0.5), vol, aon]), header=" vlow vhigh vmean scale bias biasErr volume a")
