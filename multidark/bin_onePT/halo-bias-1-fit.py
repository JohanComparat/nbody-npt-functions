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
# Tinker 05 :
aa = 1/2**0.5
bb = 0.35
cc = 0.8
deltac = 1.68
bnu = lambda nu : 1 + aa**(-0.5)/(deltac) * ( aa**0.5*(aa*nu**2.) + aa*0.5 * bb * (aa * nu**2.)**(1-cc) - (aa*nu**2.)**cc / ((aa*nu**2.)**cc + bb*(1-cc)*(1-cc/2.) ) )

bvm = lambda vm, vcut : bnu(vm/vcut)

fun = lambda vm, a, b, c :  b* (vm/c) +  a*n.e**(vm**0.5/c**0.5)
xs = n.arange(60, 2000, 1)

zmin=-0.1
zmax=0.1

sig_low, sig_high, sig_mean, scale, bias, biasErr, vol, aon,logmps = n.loadtxt(join(os.environ['MVIR_DIR'],  "halo-bias-measurement-summary.data"), unpack=True)

############
# Fit at redshift 0
############

sel0 = (aon==1)
res, cov = curve_fit(fun, sig_mean[sel0], bias[sel0], p0=(0.5, -0.5, 200), sigma=biasErr[sel0], maxfev=10000000)
"""
print res
print cov.diagonal()**0.5 

p.figure(0,(5,5))
p.axes([0.18,0.18,0.75,0.75])
p.plot(xs, fun(xs, res[0], res[1], res[2]), 'k--',label='fit')
sel = (vol==400**3.)&(aon==1)
p.errorbar(vmean[sel], bias[sel], yerr = biasErr[sel], c='k', label='SMD', fmt='none')
sel = (vol==1000**3.)&(aon==1)
p.errorbar(vmean[sel], bias[sel], yerr = biasErr[sel], c='b', label='MDPL', fmt='none')
sel = (vol==2500**3.)&(aon==1)
p.errorbar(vmean[sel], bias[sel], yerr = biasErr[sel], c='r', label='BigMD', fmt='none')
sel = (vol==4000**3.)&(aon==1)
p.errorbar(vmean[sel], bias[sel], yerr = biasErr[sel], c='m', label='HMD', fmt='none')
#-cb = p.colorbar()
#cb.set_label('z')
p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel(r'bias')
p.xlim((50, 1500))
p.ylim((0.6,4.5))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","halo-bias-z0.png"))
p.clf()

p.figure(0,(5,5))
p.axes([0.18,0.18,0.75,0.75])
sel = (vol==400**3.)&(aon==1)
p.errorbar(vmean[sel], bias[sel]/fun(vmean[sel], res[0], res[1], res[2]), yerr = biasErr[sel]/bias[sel], c='k', label='SMD', fmt='none')
sel = (vol==1000**3.)&(aon==1)
p.errorbar(vmean[sel], bias[sel]/fun(vmean[sel], res[0], res[1], res[2]), yerr = biasErr[sel]/bias[sel], c='b', label='MDPl', fmt='none')
sel = (vol==2500**3.)&(aon==1)
p.errorbar(vmean[sel], bias[sel]/fun(vmean[sel], res[0], res[1], res[2]), yerr = biasErr[sel]/bias[sel], c='r', label='BigMD', fmt='none')
sel = (vol==4000**3.)&(aon==1)
p.errorbar(vmean[sel], bias[sel]/fun(vmean[sel], res[0], res[1], res[2]), yerr = biasErr[sel]/bias[sel], c='m', label='HMD', fmt='none')
p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel(r'residual = data / model')
p.xlim((50, 1500))
p.ylim((0.8,1.2))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","halo-bias-residual-z0.png"))
p.clf()
"""

##############
# Fit at redshift trend
##############
a0, b0, c0 = res

xs = n.arange(60, 2000, 1)

a_pr = lambda zz, a1, a2 : a0 *(1+a1*zz+a2*zz**2.) # + a2 *lz**2. + a3 *lz**3. + a4 *lz**4.
b_pr = lambda zz, b1, b2 : b0 *(1+b1*zz+b2*zz**2.)# +b3*zz**3.) # + b2 *lz**2.+ b3 *lz**3.)
c_pr = lambda zz, c1, c2 : c0 *(1+c1*zz+c2*zz**2) # + b2 *lz**2.+ b3 *lz**3.)
vfG = lambda vm, zz, ps  : a_pr(zz,ps[0], ps[1])*n.e**(vm**0.5/c_pr(zz,ps[4], ps[5])**0.5) +  (b_pr(zz,ps[2], ps[3]))* (vm /c_pr(zz,ps[4], ps[5]))
#vfG = lambda vm, zz, ps  : a_pr(zz,ps[0])*n.e**(vm**0.5/c_pr(zz,ps[3], ps[4])**0.5) +  (b_pr(zz,ps[1], ps[2]))* (vm /c_pr(zz,ps[3], ps[4]))
p1=[0,0,0,0,0,0] # [0.28, -0.04, 0.7, -0.47, -0., -0.55, 0.28]
chi2fun = lambda ps : n.sum( (vfG(vmean,1/aon-1,ps) - bias)**2. / (2*biasErr)**2. )/(len(bias) - len(p1))

res2 = minimize(chi2fun, p1, method='Powell',options={'xtol': 1e-6, 'disp': True, 'maxiter' : 50000000000000})

pOpt = res2.x
cov = res2.direc
chi2perpoint = lambda ps : (vfG(vmean,1/aon-1,ps) - bias)**2. / (2*biasErr)**2. 
chi2pp = chi2perpoint(pOpt)

n.savetxt(join("..","clustering","bias0-best_fit_params.txt"),n.transpose([pOpt,cov.diagonal()**0.5]))
vs = n.arange(60,1500,2)
X,Y = n.meshgrid(vs,n.arange(zmin, zmax+0.025,0.025))
Z = vfG(X,Y,pOpt)

n.savetxt(join("..","clustering","bias0-best_fit.txt"),n.transpose([n.hstack((X)), n.hstack((Y)), n.hstack((Z))]) )

###############################
#IMPLEMENT SINGLE REDSHIFT FITS TO GET THERIGHT PARAMETRIZATION
#SAVE DEPENDECE WITH REDSHIFT POINTS B-y WRITING THEM OUT
#WRITE REDSHIFT DEPENDENCE EQUATION IN THE PAPER

#######################################################
# now plots the results of the fit
print "now plots the results of the fit"

vmax_mod, z_mod, n_mod = n.loadtxt(join("..","clustering","bias0-best_fit.txt"), unpack=True)

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(vmax_mod, n_mod, c=z_mod,s=5, marker='o',label="model", rasterized=True, vmin = zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel('bias') # slog$_{10}[ n(>M)]')
p.xlim((50, 1500))
p.ylim((0.6,4.5))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","halo-bias-zAll-model.png"))
p.clf()

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.errorbar(vmean, bias, yerr = 2*biasErr, fmt='none',elinewidth=0.5, mfc='none',ecolor='k',rasterized = True)
sc1=p.scatter(vmean, bias, c=1/aon-1,s=5, marker='o',label="model", rasterized=True, vmin = zmin, vmax = zmax)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel('bias') # slog$_{10}[ n(>M)]')
p.xlim((50, 1500))
p.ylim((0.6,4.5))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","halo-bias-zAll-data.png"))
p.clf()

print "ndof=",len(bias)-len(pOpt)
print "ndata",len(bias)
print "maxchi2distance",n.max(chi2pp)
print "Noutliers=",len((chi2pp>1.5).nonzero()[0])

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
sc1=p.scatter(vmean, 1/aon-1, c=chi2pp,s=5, marker='o', rasterized=True, vmin = 0, vmax = 1.2)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label(r"(data-model)$^2$/(err data)$^2$")
p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel(r'z') # log$_{10}[ n(>M)]')
p.xlim((50, 1500))
p.grid()
p.savefig(join("..","clustering","halo-bias-zAll-chi2pp.png"))
p.clf()




sys.exit()
p.figure(0,(5,5))
p.axes([0.18,0.18,0.75,0.75])
p.scatter(vmean, bias**0.5, c=1./aon-1, s=20, edgecolors='none')
p.plot(xs, fun(xs, res[0], res[1], res[2]), 'k--')
#p.plot(xs, bvm(xs, res2[0]]), 'r--')
cb = p.colorbar()
cb.set_label('z')
p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel(r'bias')
p.xlim((50, 3000))
p.ylim((0.6,4.5))
#p.yscale('log')
p.xscale('log')
#gl = p.legend(loc=2,fontsize=10)
#gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","ahalo-bias.png"))
p.clf()

p.figure(0,(5,5))
p.axes([0.18,0.18,0.75,0.75])
p.plot(vmean[aon==1], bias[aon==1]**0.5/fun(vmean[aon==1], res[0], res[1], res[2]), 'bo')
#p.plot(vmean[aon==1], bias[aon==1]**0.5/bvm(vmean[aon==1], res2[0]), 'r+')
p.plot(vmean[aon==1],  1+ biasErr[aon==1]*bias[aon==1]**(-0.5), 'k+')
p.plot(vmean[aon==1],  1- biasErr[aon==1]*bias[aon==1]**(-0.5), 'k+')
p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel(r'residual = data / model')
p.xlim((50, 3000))
p.ylim((0.7,1.3))
#p.yscale('log')
p.xscale('log')
#gl = p.legend(loc=2,fontsize=10)
#gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","ahalo-bias-residual.png"))
p.clf()


print res
print cov.diagonal()**0.5 / res

p.figure(0,(5,5))
p.axes([0.18,0.18,0.75,0.75])
sel = (vol==400**3.)
# print len(bias[sel])
p.errorbar(vmean[sel], bias[sel]/n.polyval(ps, vmean[sel]), yerr = biasErr[sel]/ bias[sel], c='k', label='SMD')
sel = (vol==1000**3.)
# print len(bias[sel])
p.errorbar(vmean[sel],  bias[sel]/n.polyval(ps, vmean[sel]), yerr = biasErr[sel]/ bias[sel], c='b', label='MDPl')
sel = (vol==2500**3.)
# print len(bias[sel])
p.errorbar(vmean[sel],  bias[sel]/n.polyval(ps, vmean[sel]), yerr = biasErr[sel]/ bias[sel], c='r', label='BigMD')
sel = (vol==4000**3.)
# print len(bias[sel])
p.errorbar(vmean[sel],  bias[sel]/n.polyval(ps, vmean[sel]), yerr = biasErr[sel]/ bias[sel], c='m', label='HMD')

p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel(r'$b^2(V_{max})$/model')
p.xlim((50, 3000))
p.ylim((0.9,1.1))
#p.yscale('log')
p.xscale('log')
gl = p.legend(loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","a_1.00000","ahalo-bias-z0-residual-zoom.png"))
p.clf()


p.figure(0,(5,5))
p.axes([0.18,0.18,0.75,0.75])
sel = (vol==400**3.)
# print len(bias[sel])
p.errorbar(vmean[sel], bias[sel], yerr = biasErr[sel], c='k', label='SMD')
sel = (vol==1000**3.)
# print len(bias[sel])
p.errorbar(vmean[sel], bias[sel], yerr = biasErr[sel], c='b', label='MDPl')
sel = (vol==2500**3.)
# print len(bias[sel])
p.errorbar(vmean[sel], bias[sel], yerr = biasErr[sel], c='r', label='BigMD')
sel = (vol==4000**3.)
# print len(bias[sel])
p.errorbar(vmean[sel], bias[sel], yerr = biasErr[sel], c='m', label='HMD')
p.plot(vmean, n.polyval(ps, vmean), 'r--', lw=2,label='fit')

p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel(r'$b^2(V_{max}$)')
p.xlim((50, 3000))
#p.ylim((0.1,100))
p.yscale('log')
p.xscale('log')
gl = p.legend(loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","a_1.00000","ahalo-bias-z0.png"))
p.clf()

	
	
sys.exit()

ps = n.polyfit(vmean**2., bias, degree)#, w = 1./(biasErr))
n.savetxt("fit-halo-bias2-vmax2.data",ps)

p.figure(0,(5,5))
p.axes([0.18,0.18,0.75,0.75])
sel = (vol==400**3.)
# print len(bias[sel])
p.errorbar(vmean[sel]**2., bias[sel]/n.polyval(ps, vmean[sel]**2.), yerr = biasErr[sel]/ bias[sel], c='k', label='SMD')
sel = (vol==1000**3.)
# print len(bias[sel])
p.errorbar(vmean[sel]**2.,  bias[sel]/n.polyval(ps, vmean[sel]**2.), yerr = biasErr[sel]/ bias[sel], c='b', label='MDPl')
sel = (vol==2500**3.)
# print len(bias[sel])
p.errorbar(vmean[sel]**2.,  bias[sel]/n.polyval(ps, vmean[sel]**2.), yerr = biasErr[sel]/ bias[sel], c='r', label='BigMD')
sel = (vol==4000**3.)
# print len(bias[sel])
p.errorbar(vmean[sel]**2.,  bias[sel]/n.polyval(ps, vmean[sel]**2.), yerr = biasErr[sel]/ bias[sel], c='m', label='HMD')

p.xlabel(r'$V_{max}$ (km/s)')
p.ylabel(r'$b^2(V_{max})$/model')
#p.xlim((50, 3000))
p.ylim((0.7,1.3))
#p.yscale('log')
p.xscale('log')
gl = p.legend(loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","a_1.00000","bhalo-bias-vm2-z0-residual.png"))
p.clf()


p.figure(0,(5,5))
p.axes([0.18,0.18,0.75,0.75])
sel = (vol==400**3.)
# print len(bias[sel])
p.errorbar(vmean[sel]**2., bias[sel], yerr = biasErr[sel], c='k', label='SMD')
sel = (vol==1000**3.)
# print len(bias[sel])
p.errorbar(vmean[sel]**2., bias[sel], yerr = biasErr[sel], c='b', label='MDPl')
sel = (vol==2500**3.)
# print len(bias[sel])
p.errorbar(vmean[sel]**2., bias[sel], yerr = biasErr[sel], c='r', label='BigMD')
sel = (vol==4000**3.)
# print len(bias[sel])
p.errorbar(vmean[sel]**2., bias[sel], yerr = biasErr[sel], c='m', label='HMD')
p.plot(vmean**2., n.polyval(ps, vmean**2.), 'r--', lw=2,label='fit')

p.xlabel(r'$V^2_{max}$ (km/s)')
p.ylabel(r'$b^2(V_{max}$)')
#p.xlim((50, 3000))
#p.ylim((0.1,100))
p.yscale('log')
p.xscale('log')
gl = p.legend(loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig(join("..","clustering","a_1.00000","bhalo-bias-vm2-z0.png"))
p.clf()

sys.exit()



































for vmin in vmins:
	list44 = n.array(glob.glob(join("..","MD_1Gpc","halo_bias","clustering","hlist_1.00000_vmax_"+vmin+"*rmax_050_xiR.pkl")))
	list42 = n.array(glob.glob(join("..","MD_0.4Gpc","halo_bias","clustering","hlist_1.00000_vmax_"+vmin+"*rmax_050_xiR.pkl")))
	list41 = n.array(glob.glob(join("..","MD_2.5Gpc","halo_bias","clustering","hlist_80_vmax_"+vmin+"*rmax_140_xiR.pkl")))
	list43 = n.array(glob.glob(join("..","MD_4Gpc","halo_bias","clustering","hlist_128_vmax_"+vmin+"*rmax_140_xiR.pkl")))
	list40=n.hstack((list41, list42, list43, list44))
	list40.sort()
	# print list40
	p.figure(0,(11,6))
	p.axes([0.15,0.15,0.6,0.75])
	for ii in range(len(list40)):
		f=open(list40[ii],'r')
		bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbinsL, vbinsH = cPickle.load(f)
		f.close()
		if list40[ii].split('\\')[1] == "MD_0.4Gpc":
			color = 'w'
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
	p.ylim((-50,100))
	p.title(str(n.round(vbinsL))+"<vmax<"+str(n.round(vbinsH))+" z=0")
	#p.yscale('log')
	#p.xscale('log')
	gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
	gl.set_frame_on(False)
	p.grid()
	p.savefig(join("..","clustering","a_1.00000","xi-MD-"+vmin+".png"))
	p.clf()


for vmin in vmins:
	list44 = n.array(glob.glob(join("..","MD_1Gpc","halo_bias","clustering","hlist_1.00000_vmax_"+vmin+"*rmax_015_xiR.pkl")))
	list42 = n.array(glob.glob(join("..","MD_0.4Gpc","halo_bias","clustering","hlist_1.00000_vmax_"+vmin+"*rmax_015_xiR.pkl")))
	list41 = n.array(glob.glob(join("..","MD_2.5Gpc","halo_bias","clustering","hlist_80_vmax_"+vmin+"*rmax_015_xiR.pkl")))
	list43 = n.array(glob.glob(join("..","MD_4Gpc","halo_bias","clustering","hlist_128_vmax_"+vmin+"*rmax_015_xiR.pkl")))
	list40=n.hstack((list41, list42, list43, list44))
	list40.sort()
	# print list40
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
# print list40

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
# print list40

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
# print list04

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
