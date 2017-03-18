import astropy.units as uu
import astropy.cosmology as co
aa = co.Planck13
#aah = co.FlatLambdaCDM(H0=100.0 *uu.km / (uu.Mpc *uu.s), Om0=0.307, Tcmb0=2.725 *uu.K, Neff=3.05, m_nu=[ 0.  ,  0. ,   0.06]*uu.eV, Ob0=0.0483)
#rhom = aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aah.critical_density0.to(uu.solMass*uu.Mpc**-3).value

limits_04 = [70, 3000]
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


##########################################################3
##########################################################3
##########################################################3
# Z=1
##########################################################3
##########################################################3
##########################################################3



list04 = n.array(["D:\data\MultiDark\MD_0.4Gpc\halo_bias\clustering\hlist_0.50000_vmax_190.0_213.18_xiR.pkl",
"D:\data\MultiDark\MD_0.4Gpc\halo_bias\clustering\hlist_0.50000_vmax_213.18_239.2_xiR.pkl",
"D:\data\MultiDark\MD_0.4Gpc\halo_bias\clustering\hlist_0.50000_vmax_239.2_268.38_xiR.pkl",
"D:\data\MultiDark\MD_0.4Gpc\halo_bias\clustering\hlist_0.50000_vmax_268.38_301.13_xiR.pkl"])

list25 = n.array(["D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_190.0_213.18_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_213.18_239.2_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_239.2_268.38_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_268.38_301.13_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_301.13_337.87_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_337.87_379.1_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_379.1_425.36_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_425.36_477.26_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_477.26_535.49_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_535.49_600.83_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_600.83_674.15_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_674.15_756.4_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_756.4_848.7_xiR.pkl",
"D:\data\MultiDark\MD_2.5Gpc\halo_bias\clustering\hlist_11_vmax_848.7_952.26_xiR.pkl"])

list10 = n.array(["D:\data\MultiDark\MD_1Gpc\halo_bias\clustering\hlist_0.49220_vmax_190.0_213.18_xiR.pkl",
"D:\data\MultiDark\MD_1Gpc\halo_bias\clustering\hlist_0.49220_vmax_213.18_239.2_xiR.pkl",
"D:\data\MultiDark\MD_1Gpc\halo_bias\clustering\hlist_0.49220_vmax_239.2_268.38_xiR.pkl",
"D:\data\MultiDark\MD_1Gpc\halo_bias\clustering\hlist_0.49220_vmax_268.38_301.13_xiR.pkl",
"D:\data\MultiDark\MD_1Gpc\halo_bias\clustering\hlist_0.49220_vmax_301.13_337.87_xiR.pkl",
"D:\data\MultiDark\MD_1Gpc\halo_bias\clustering\hlist_0.49220_vmax_337.87_379.1_xiR.pkl",
"D:\data\MultiDark\MD_1Gpc\halo_bias\clustering\hlist_0.49220_vmax_379.1_425.36_xiR.pkl",
"D:\data\MultiDark\MD_1Gpc\halo_bias\clustering\hlist_0.49220_vmax_425.36_477.26_xiR.pkl",
"D:\data\MultiDark\MD_1Gpc\halo_bias\clustering\hlist_0.49220_vmax_477.26_535.49_xiR.pkl"])

list04.sort()
list10.sort()
list25.sort()

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
p.plot(xrr,20*xrr**(-1.8),'k--',lw=2)
p.axvline(3)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel('xi(MDPL) (r)')
p.xlim((0.5,20))
p.ylim((0.01,100))
p.yscale('log')
p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig('xi-SMDPL-z1.png')
p.clf()

p.figure(0,(11,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list10)): #[::2]:
	f=open(list10[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	fv_04 = 1000**3. / volume_04 
	xi04V = fv_04*(xis_04+1)-1.
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	p.plot(rr,xi04V,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))


xrr = n.arange(0,50,0.5)
p.plot(xrr,20*xrr**(-1.8),'k--',lw=2)
p.axvline(3)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel('xi(MDPL) (r)')
p.xlim((0.5,20))
p.ylim((0.01,100))
p.yscale('log')
p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig('xi-MDPL-z1.png')
p.clf()


p.figure(0,(11,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list25)): #[::2]:
	f=open(list25[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	fv_04 = 2500**3. / volume_04 
	xi04V = fv_04*(xis_04+1)-1.
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	p.plot(rr,xi04V,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))

xrr = n.arange(0,50,0.5)
p.plot(xrr,20*xrr**(-1.8),'k--',lw=2)
p.axvline(3)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel('xi(BigMD) (r)')
p.xlim((0.5,20))
p.ylim((0.01,100))
p.yscale('log')
p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig('xi-BigMD-z1.png')
p.clf()




p.figure(0,(10,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list04)):
	f=open(list04[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	#p.plot(rr,xis,label="0.4 190-213")
	f=open(list25[ii],'r')
	bin_xi3D_25,xis_25, DR_25, volume_25, dV_25, pairCount_25, pairs_25, Ntotal_25, nD_25, nR_25, vbinsL_25, vbinsH_25 = cPickle.load(f)
	f.close()
	ok = (rr>3)&(rr<7)
	fv_04 = 400.**3. / volume_04 
	fv_25 = 2500.**3. / volume_25
	xi25V = fv_25*(xis_25[ok]+1)-1.
	xi04V = fv_04*(xis_04[ok]+1)-1.
	y = xi25V / xi04V
	print y
	p.plot(rr[ok],y,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))


p.axvline(3)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel(' xi(BigMD)/xi(SMDPL) (r)')
p.xlim((2.9,7.5))
p.ylim((0.85,1.15))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig('xi-ratio-BigMD-SMDPL-z1.png')
p.show()


p.figure(0,(10,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list04)):
	f=open(list04[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	#p.plot(rr,xis,label="0.4 190-213")
	f=open(list10[ii],'r')
	bin_xi3D_25,xis_25, DR_25, volume_25, dV_25, pairCount_25, pairs_25, Ntotal_25, nD_25, nR_25, vbinsL_25, vbinsH_25 = cPickle.load(f)
	f.close()
	ok = (rr>3)&(rr<7)
	fv_04 = 400.**3. / volume_04 
	fv_25 = 1000.**3. / volume_25
	xi25V = fv_25*(xis_25[ok]+1)-1.
	xi04V = fv_04*(xis_04[ok]+1)-1.
	y = xi25V / xi04V
	print y
	p.plot(rr[ok],y,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))


p.axvline(3)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel(' xi(MDPL)/xi(SMDPL) (r)')
p.xlim((2.9,7.5))
p.ylim((0.85,1.15))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig('xi-ratio-MDPL-SMDPL-z1.png')
p.show()


p.figure(0,(10,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list10)):
	f=open(list10[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	#p.plot(rr,xis,label="0.4 190-213")
	f=open(list25[ii],'r')
	bin_xi3D_25,xis_25, DR_25, volume_25, dV_25, pairCount_25, pairs_25, Ntotal_25, nD_25, nR_25, vbinsL_25, vbinsH_25 = cPickle.load(f)
	f.close()
	ok = (rr>3)&(rr<7)
	fv_04 = 1000.**3. / volume_04 
	fv_25 = 2500.**3. / volume_25
	xi25V = fv_25*(xis_25[ok]+1)-1.
	xi04V = fv_04*(xis_04[ok]+1)-1.
	y = xi25V / xi04V
	print y
	p.plot(rr[ok],y,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))


p.axvline(3)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel(' xi(BigMD)/xi(MDPL) (r)')
p.xlim((2.9,7.5))
p.ylim((0.85,1.15))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig('xi-ratio-BigMD-MDPL-z1.png')
p.show()


sys.exit()



##########################################################3
##########################################################3
##########################################################3
# Z=1.5
##########################################################3
##########################################################3
##########################################################3












list10 = n.array(["MD_1Gpc\halo_bias\clustering\hlist_0.40320_vmax_190.0_213.18_xiR.pkl"
,"MD_1Gpc\halo_bias\clustering\hlist_0.40320_vmax_213.18_239.2_xiR.pkl"
,"MD_1Gpc\halo_bias\clustering\hlist_0.40320_vmax_239.2_268.38_xiR.pkl"
,"MD_1Gpc\halo_bias\clustering\hlist_0.40320_vmax_268.38_301.13_xiR.pkl"
,"MD_1Gpc\halo_bias\clustering\hlist_0.40320_vmax_301.13_337.87_xiR.pkl"
,"MD_1Gpc\halo_bias\clustering\hlist_0.40320_vmax_337.87_379.1_xiR.pkl"
,"MD_1Gpc\halo_bias\clustering\hlist_0.40320_vmax_379.1_425.36_xiR.pkl"
,"MD_1Gpc\halo_bias\clustering\hlist_0.40320_vmax_425.36_477.26_xiR.pkl"
,"MD_1Gpc\halo_bias\clustering\hlist_0.40320_vmax_477.26_535.49_xiR.pkl"])

list25 = n.array(["MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_190.0_213.18_xiR.pkl"
,"MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_213.18_239.2_xiR.pkl"
,"MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_239.2_268.38_xiR.pkl"
,"MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_268.38_301.13_xiR.pkl"
,"MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_301.13_337.87_xiR.pkl"
,"MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_337.87_379.1_xiR.pkl"
,"MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_379.1_425.36_xiR.pkl"
,"MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_425.36_477.26_xiR.pkl"
,"MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_477.26_535.49_xiR.pkl"])

list25All = n.array(glob.glob("MD_2.5Gpc\halo_bias\clustering\hlist_10_vmax_*_*_xiR.pkl"))

list10.sort()
list25.sort()

list25All.sort()

f=open(list10[-1],'r')
bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
f.close()
rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.

#p.plot(rr,xis,label="0.4 190-213")
f=open(list25[-1],'r')
bin_xi3D_25,xis_25, DR_25, volume_25, dV_25, pairCount_25, pairs_25, Ntotal_25, nD_25, nR_25, vbinsL_25, vbinsH_25 = cPickle.load(f)
f.close()

ok = (rr>2)&(rr<7)
fv_04 = 1000**3. / volume_04 
fv_25 = 2500.**3. / volume_25
xi25V = fv_25*(xis_25[ok]+1)-1.
xi04V = fv_04*(xis_04[ok]+1)-1.
xi25V / xi04V

fc_04 = pairCount_04 / (nR_04*nD_04 -nD_04 *(nD_04 -1)/2.) 
fc_25 = pairCount_25/ (nR_25*nD_25 -nD_25 *(nD_25 -1)/2.) 
xi25C = fc_25*(xis_25[ok]+1)-1.
xi04C = fc_04*(xis_04[ok]+1)-1.
xi25C / xi04C



y= (DR_25* volume_04)/ (DR_04 * volume_25)
y= (DR_25* dV_04)/ (DR_04 * dV_25)
print y
p.plot(rr,y,label="Npairs(2.5)/(0.4) 190<vmax<213")




p.figure(0,(11,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list10)): #[::2]:
	f=open(list10[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	fv_04 = 1000**3. / volume_04 
	xi04V = fv_04*(xis_04+1)-1.
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	p.plot(rr,xi04V,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))


xrr = n.arange(0,50,0.5)
p.plot(xrr,60*xrr**(-1.8),'k--',lw=2)
p.axvline(2)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel('xi(MDPL) (r)')
p.xlim((0.5,20))
p.ylim((0.01,100))
p.yscale('log')
p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig('xi-MDPL.png')
p.clf()


p.figure(0,(11,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list25All)): #[::2]:
	f=open(list25All[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	fv_04 = 2500**3. / volume_04 
	xi04V = fv_04*(xis_04+1)-1.
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	p.plot(rr,xi04V,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))

xrr = n.arange(0,50,0.5)
p.plot(xrr,60*xrr**(-1.8),'k--',lw=2)
p.axvline(2)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel('xi(BigMD) (r)')
p.xlim((0.5,20))
p.ylim((0.01,100))
p.yscale('log')
p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig('xi-BigMD.png')
p.clf()

p.figure(0,(10,6))
p.axes([0.17,0.17,0.6,0.8])

for ii in range(len(list25))[:6]:
	f=open(list10[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	#p.plot(rr,xis,label="0.4 190-213")
	f=open(list25[ii],'r')
	bin_xi3D_25,xis_25, DR_25, volume_25, dV_25, pairCount_25, pairs_25, Ntotal_25, nD_25, nR_25, vbinsL_25, vbinsH_25 = cPickle.load(f)
	f.close()
	ok = (rr>2)&(rr<7)
	fv_04 = 1000.**3. / volume_04 
	fv_25 = 2500.**3. / volume_25
	xi25V = fv_25*(xis_25[ok]+1)-1.
	xi04V = fv_04*(xis_04[ok]+1)-1.
	y = xi25V / xi04V
	print y
	p.plot(rr[ok],y,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)))


for ii in range(len(list25))[6:-1]:
	f=open(list10[ii],'r')
	bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
	f.close()
	rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
	#p.plot(rr,xis,label="0.4 190-213")
	f=open(list25[ii],'r')
	bin_xi3D_25,xis_25, DR_25, volume_25, dV_25, pairCount_25, pairs_25, Ntotal_25, nD_25, nR_25, vbinsL_25, vbinsH_25 = cPickle.load(f)
	f.close()
	ok = (rr>2)&(rr<7)
	fv_04 = 1000.**3. / volume_04 
	fv_25 = 2500.**3. / volume_25
	xi25V = fv_25*(xis_25[ok]+1)-1.
	xi04V = fv_04*(xis_04[ok]+1)-1.
	y = xi25V / xi04V
	print y
	p.plot(rr[ok],y,label= str(n.round(vbinsL_04))+"<vmax<"+str(n.round(vbinsH_04)),lw=2,ls='dashed')


p.axvline(2)
p.axvline(7)
p.xlabel('r Mpc/h')
p.ylabel(' xi(BigMD)/xi(MDPL) (r)')
p.xlim((1.9,7.5))
p.ylim((0.9,1.1))
#p.yscale('log')
#p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.savefig('xi-ratio-BigMD-MDPL.png')
p.show()


sys.exit()

p.figure(0,(10,6))
p.axes([0.17,0.17,0.6,0.8])
f=open("MD_0.4Gpc\halo_bias\clustering\hlist_1.00000_vmax_190.0_213.18_xiR.pkl",'r')
bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
f.close()
rr = (bin_xi3D_04[1:] + bin_xi3D_04[:-1])/2.
#p.plot(rr,xis,label="0.4 190-213")
f=open("MD_2.5Gpc\halo_bias\clustering\hlist_80_vmax_190.0_213.18_xiR.pkl",'r')
bin_xi3D_25,xis_25, DR_25, volume_25, dV_25, pairCount_25, pairs_25, Ntotal_25, nD_25, nR_25, vbinsL_25, vbinsH_25 = cPickle.load(f)
f.close()
y= (DR_25* volume_04)/ (DR_04 * volume_25)
y= (DR_25* dV_04)/ (DR_04 * dV_25)
print y
p.plot(rr,y,label="Npairs(2.5)/(0.4) 190<vmax<213")

f=open("MD_0.4Gpc\halo_bias\clustering\hlist_1.00000_vmax_213.18_239.2_xiR.pkl",'r')
bin_xi3D_04,xis_04, DR_04, volume_04, dV_04, pairCount_04, pairs_04, Ntotal_04, nD_04, nR_04, vbinsL_04, vbinsH_04 = cPickle.load(f)
f.close()
f=open("MD_2.5Gpc\halo_bias\clustering\hlist_80_vmax_213.18_239.2_xiR.pkl",'r')
bin_xi3D_25,xis_25, DR_25, volume_25, dV_25, pairCount_25, pairs_25, Ntotal_25, nD_25, nR_25, vbinsL_25, vbinsH_25= cPickle.load(f)
f.close()
y=(DR_25* volume_04)/ (DR_04 * volume_25)
y= (DR_25* dV_04)/ (DR_04 * dV_25)
print y
p.plot(rr,y,label="Npairs(2.5)/(0.4) 213<vmax<240")

p.xlabel('r Mpc/h')
p.ylabel(' r xi(2.5Gpc)/xi(0.4Gpc)(r)')
#p.yscale('log')
#p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.show()


sys.exit()


ll = n.array( glob.glob(join(dir_10,"halo_bias","clustering","*.pkl")))#dir_25_b,"*.pkl" ) ))
ll.sort()
scaleFact = n.empty(len(ll))
vmin = n.empty(len(ll))
vmax = n.empty(len(ll))
zz = n.empty(len(ll))
for ii in range(len(ll)):
	scaleFact[ii] = float(ll[ii].split('\\')[-1].split('_')[1])
	vmin[ii] = float(ll[ii].split('\\')[-1].split('_')[3])
	vmax[ii] = float(ll[ii].split('\\')[-1].split('_')[4])
	zz[ii] = 1./scaleFact[ii] -1.

f=open(,"MD_2.5Gpc\halo_bias\clustering\hlist_46_vmax_200.0_224.4_xiR.pkl",'r')
bin_xi3D_25, xi_25_200 = cPickle.load(f)
f.close()
f=open(,"MD_2.5Gpc\halo_bias\clustering\hlist_46_vmax_224.4_251.79_xiR.pkl",'r')
bin_xi3D_25, xi_25_224 = cPickle.load(f)
f.close()
rr_25 = (bin_xi3D_25[1:] + bin_xi3D_25[:-1])/2.


sel = (zz>2.9)
xi = n.empty((len(ll[sel]),20))
for ii, file in enumerate(ll[sel]):
	f=open(file,'r')
	bin_xi3D, xi[ii] = cPickle.load(f)
	f.close()

rr = (bin_xi3D[1:] + bin_xi3D[:-1])/2.
	
p.figure(0,(10,6))
p.axes([0.17,0.17,0.6,0.8])
p.plot(rr_25,xi_25_200,ls='dashed', lw=2,label="L2.5 z0.3 vmax200 ")
#p.plot(rr_25,xi_25_224,ls='dashed', lw=2,label="L2.5 z0.3 vmax224 ")
for ii in n.arange(len(ll[sel])):
	p.plot(rr,xi[ii],label=str(n.round(vmin[sel][ii])))

p.xlabel('r Mpc/h')
p.ylabel(' r xi(r)')
p.yscale('log')
p.xscale('log')
gl = p.legend(bbox_to_anchor=(1.05, 1), loc=2,fontsize=10)
gl.set_frame_on(False)
p.grid()
p.show()
	