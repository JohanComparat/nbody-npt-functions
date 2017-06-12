import astropy.units as uu
import astropy.cosmology as co
aa = co.Planck13
from astropy.io.fits import fits
#aah = co.FlatLambdaCDM(H0=100.0 *uu.km / (uu.Mpc *uu.s), Om0=0.307, Tcmb0=2.725 *uu.K, Neff=3.05, m_nu=[ 0.  ,  0. ,   0.06]*uu.eV, Ob0=0.0483)
#rhom = aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aah.critical_density0.to(uu.solMass*uu.Mpc**-3).value

limits_04 = [100, 3000]
limits_10 = [200, 3000]
limits_25 = [400, 3000]
limits_40 = [650, 5000]

from mpl_toolkits.mplot3d import Axes3D
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
dir_04 = join(dir,"MD_0.4Gpc")
dir_10 = join(dir,"MD_1Gpc")
dir_25 = join(dir,"MD_2.5Gpc")
dir_40 = join(dir,"MD_4Gpc")

dirDF = join(dir_10, "density_field")

mocks = n.array(glob.glob(join(dirDF,"Box*.gz")))

jj=0
datLRG = fits.open(mocks[jj])[1].data
jj=1
datQSO = fits.open(mocks[jj])[1].data

p.plot(datLRG['DF_N1'], datLRG['Vmax'],'r+')
p.plot(datQSO['DF_N1'], datQSO['Vmax'],'b+')
p.xlabel(r'$\delta$')
p.ylabel(r'$V_{max}$')
p.xscale('log')
p.yscale('log')
p.show()

v = datLRG['vx']**2.+datLRG['vy']**2.+datLRG['vz']**2.
p.plot(v**0.5, datLRG['DF'],'r+')
v = datQSO['vx']**2.+datQSO['vy']**2.+datQSO['vz']**2.
p.plot(v**0.5, datQSO['DF'],'b+')
p.ylabel(r'$\delta$')
p.xlabel(r'$V$')
p.xscale('log')
p.yscale('log')
p.show()






zList_file =  join(dir_25, "redshift_list.txt") 
n0,z0,a0 = n.loadtxt(zList_file,unpack=True)
n1 = n0.astype(int)+1
n2 = n1.astype(str)
snap2Z = {n2[i].zfill(3): z0[i] for i in range(len(n2))}

j=79
Rs, xiR  = n.loadtxt(join(dir,"Pk_DM_CLASS","multidark00_z"+n2[j]+"_xi.dat"), unpack = True)


vmins = ["65.0", "72.9", "81.8", "91.8", "103.0", "115.5", "129.6", "145.5", "163.2", "183.1", "205.5", "230.6","258.7", "290.3", "325", "365.5", "410.1", "460.1", "516.3", "579.3", "650.0", "729.3", "818.3", "918.1", "1030.1", "1155.8", "1296.9"]

for vmin in vmins:
	list42 = n.array(glob.glob(join("..","MD_*Gpc","halo_bias","clustering","hlist_1.00000_vmax_"+vmin+"*rmax_140_xiR.pkl")))
	list41 = n.array(glob.glob(join("..","MD_2.5Gpc","halo_bias","clustering","hlist_80_vmax_"+vmin+"*rmax_140_xiR.pkl")))
	list43 = n.array(glob.glob(join("..","MD_4Gpc","halo_bias","clustering","hlist_128_vmax_"+vmin+"*rmax_140_xiR.pkl")))
	list40=n.hstack((list41, list42, list43))
	list40.sort()
	print list40
	if len(list40)==2:
		vols=[1000. , 2500.]
	if len(list40)==3:
		vols=[400., 1000. , 2500.]

	p.figure(0,(11,6))
	p.axes([0.15,0.15,0.6,0.75])
	for ii in range(len(list40)):
		f=open(list40[ii],'r')
		bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbinsL, vbinsH = cPickle.load(f)
		f.close()
		if list40[ii].split('\\')[1] == "MD_0.4Gpc":
			color = 'k'
			volume = 400**3.
		if list40[ii].split('\\')[1] == "MD_1Gpc":
			color = 'b'
			volume = 1000**3.
		if list40[ii].split('\\')[1] == "MD_2.5Gpc":
			color = 'r'
			volume = 2500**3.
		if list40[ii].split('\\')[1] == "MD_4Gpc":
			color = 'm'
			volume = 4000**3.
		
		DR_rb = DR[::2][:-1] + DR[1::2]
		dV_rb =  dV[::2][:-1] + dV[1::2]
		xi_rb = DR_rb*volume/(dV_rb * pairCount) -1.
		rr = (bin_xi3D[1:] + bin_xi3D[:-1])/2.
		rr_rb = bin_xi3D[::2][1:]
		p.plot(rr_rb, rr_rb*rr_rb*xi_rb,label= list40[ii].split('\\')[1], c = color)
		
	p.plot(Rs,Rs*Rs*xiR,'b--',label='DM linear theory')
	p.xlabel('r Mpc/h')
	p.ylabel(r'$r^2 \xi$(MD) (r)')
	p.xlim((0,200))
	p.ylim((-400,400))
	p.title(str(n.round(vbinsL))+"<vmax<"+str(n.round(vbinsH))+" z=0")
	#p.yscale('log')
	#p.xscale('log')
	gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
	gl.set_frame_on(False)
	p.grid()
	p.savefig(join("..","clustering","a_1.00000","xi-MD-"+vmin+".png"))
	p.clf()

for vmin in vmins:
	list42 = n.array(glob.glob(join("..","MD_*Gpc","halo_bias","clustering","hlist_1.00000_vmax_"+vmin+"*rmax_015_xiR.pkl")))
	list41 = n.array(glob.glob(join("..","MD_2.5Gpc","halo_bias","clustering","hlist_80_vmax_"+vmin+"*rmax_015_xiR.pkl")))
	list43 = n.array(glob.glob(join("..","MD_4Gpc","halo_bias","clustering","hlist_128_vmax_"+vmin+"*rmax_015_xiR.pkl")))
	list40=n.hstack((list41, list42, list43))
	list40.sort()
	print list40
	if len(list40)==2:
		vols=[1000. , 2500.]
	if len(list40)==3:
		vols=[400., 1000. , 2500.]

	p.figure(0,(11,6))
	p.axes([0.15,0.15,0.6,0.75])
	for ii in range(len(list40)):
		f=open(list40[ii],'r')
		bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbinsL, vbinsH = cPickle.load(f)
		f.close()
		if list40[ii].split('\\')[1] == "MD_0.4Gpc":
			color = 'k'
			volume = 400**3.
		if list40[ii].split('\\')[1] == "MD_1Gpc":
			color = 'b'
			volume = 1000**3.
		if list40[ii].split('\\')[1] == "MD_2.5Gpc":
			color = 'r'
			volume = 2500**3.
		if list40[ii].split('\\')[1] == "MD_4Gpc":
			color = 'm'
			volume = 4000**3.
		
		xi = DR*volume/(dV * pairCount) -1.
		rr = (bin_xi3D[1:] + bin_xi3D[:-1])/2.
		p.plot(rr, rr*xi,label= list40[ii].split('\\')[1], c = color)
		
	p.plot(Rs,Rs*xiR,'b--',label='DM linear theory')
	p.xlabel('r Mpc/h')
	p.ylabel(r'$r \xi$(MD) (r)')
	p.xlim((0.01,20))
	p.ylim((1.,200))
	p.title(str(n.round(vbinsL))+"<vmax<"+str(n.round(vbinsH))+" z=0")
	p.yscale('log')
	p.xscale('log')
	gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
	gl.set_frame_on(False)
	p.grid()
	p.savefig(join("..","clustering","a_1.00000","xi-lt20-MD-"+vmin+".png"))
	p.clf()

sys.exit()
##########################################################3
##########################################################3
##########################################################3
# Z=1
##########################################################3
##########################################################3
##########################################################3

list40 = n.array(glob.glob(join("..","MD_1Gpc","halo_bias","clustering","hlist_1.0*_vmax_*_xiR.pkl")))

list40.sort()
print list40

p.figure(0,(11,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list40)):
	f=open(list40[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	fv_04 = 1000**3. / volume_04 
	xi04V = fv_04*(xis_04+1)-1.
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	p.plot(rr,rr*rr*xi04V,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))


p.xlabel('r Mpc/h')
p.ylabel(r'$r^2 \xi$(BigMDPL) (r)')
p.xlim((0,200))
p.ylim((-1,200))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","xi-MDPL.png"))
p.show()


list40 = n.array(glob.glob(join("..","MD_2.5Gpc","halo_bias","clustering","hlist_80*_vmax_*_xiR.pkl")))

list40.sort()
print list40

p.figure(0,(11,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list40))[:-3][::2]:
	f=open(list40[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	fv_04 = 2500**3. / volume_04 
	xi04V = fv_04*(xis_04+1)-1.
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	p.plot(rr,rr*rr*xi04V,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))


p.xlabel('r Mpc/h')
p.ylabel(r'$r^2 \xi$(BigMDPL) (r)')
p.xlim((0,200))
p.ylim((-1,200))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","xi-BigMDPL.png"))
p.show()

list04 = n.array(glob.glob(join("..","MD_0.4Gpc","halo_bias","clustering","hlist_1.00*_vmax_*_xiR.pkl")))

list04.sort()
print list04

p.figure(0,(11,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list04)): #[::2]:
	f=open(list04[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	fv_04 = 400**3. / volume_04 
	xi04V = fv_04*(xis_04+1)-1.
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	p.plot(rr,xi04V,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))


xrr = n.arange(0,50,0.5)
#p.plot(xrr,20*xrr**(-1.8),'k--',lw=2)
p.axvline(3)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel('xi(MDPL) (r)')
p.xlim((0.1,15))
p.ylim((0.1,200))
p.yscale('log')
p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","xi-SMDPL.png"))
p.show()
