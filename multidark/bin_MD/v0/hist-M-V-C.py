import pylab as p
import glob
import numpy as n
import astropy.cosmology as co
aa=co.Planck13
import astropy.units as uu
import cPickle
import sys
from scipy.interpolate import interp1d

volume=10**9

snList=[]
for i in range(100):
	snL=glob.glob("/Volumes/data/BigMD/1Gpc_3840_Planck1/MFMC/*"+str(i)+"*MF.hist.dat")
	if len(snL)==5:
		snList.append(i)

snList=n.array(snList)
snapNum,snapZ,snapA=n.loadtxt("/Volumes/data/BigMD/1Gpc_3840_Planck1/redshift-snapshot.list",unpack=True)

snL=glob.glob("/data2/DATA/eBOSS/Multidark-properties/MDPL/*.cat.gz")

numtoZ=interp1d(snapNum,snapZ)

js=[0,1,3,10,30,50]
for jj in js:
	print jj
	snList[jj]
	snL=glob.glob("/Volumes/data/BigMD/1Gpc_3840_Planck1/MFMC/*"+ str(snList[jj])+ "*MF.hist.dat")
	nnC=n.empty([6,799])*0.
	nnS=n.empty([6,799])*0.
	for iii in range(len(snL)):
		mmin,mmax,nnC[iii],nnS[iii]=n.loadtxt(snL[iii],unpack=True)

	nC=n.sum(nnC,axis=0)
	nS=n.sum(nnS,axis=0)
	mass=(mmin+mmax)/2.
	dLogM=0.01 # n.median(mmax-mmin)
	sel=(nC>=10)&(nC<3840**3)
	dmf=interp1d(mass[sel],nC[sel]/dLogM/volume)

	p.plot(mass[sel],nC[sel]/dLogM/volume,ls='None',marker='+',label='halos z='+str(numtoZ(snList[jj])))
	#p.plot(mass[sel],nS[sel]/dLogM/volume,'b+',label='sat')
	
p.grid()
p.xlabel('log halo mass')
p.ylabel('dN/dlogM/V $h^4$ Mpc$^{-3}M_\odot^{-1}$')
p.yscale('log')
p.ylim((1e-7,1))
p.xlim((10,16))
p.legend(fontsize=9)
p.savefig("1Gpc_3840_Planck1/MFplot/massfunction-evolution.pdf")
p.clf()



js=[0,1,3,10,30,50]
for jj in js:
	print jj
	snList[jj]
	snL=glob.glob("/Volumes/data/BigMD/1Gpc_3840_Planck1/MFMC/*"+ str(snList[jj])+ "*MF.hist.dat")
	nnC=n.empty([6,799])*0.
	nnS=n.empty([6,799])*0.
	for iii in range(len(snL)):
		mmin,mmax,nnC[iii],nnS[iii]=n.loadtxt(snL[iii],unpack=True)

	nC=n.sum(nnC,axis=0)
	nS=n.sum(nnS,axis=0)
	mass=(mmin+mmax)/2.
	dLogM=0.01 # n.median(mmax-mmin)
	sel=(nC>=10)&(nC<3840**3)

	p.plot(mass[sel],nC[sel]/dLogM/volume,ls='None',marker='+',label='halos z='+str(numtoZ(snList[jj])))
	#p.plot(mass[sel],nS[sel]/dLogM/volume,'b+',label='sat')
	
p.grid()
p.xlabel('log halo mass')
p.ylabel('dN/dlogM/V $h^4$ Mpc$^{-3}M_\odot^{-1}$')
p.yscale('log')
p.ylim((1e-7,1))
p.xlim((10,16))
p.legend(fontsize=9)
p.savefig("1Gpc_3840_Planck1/MFplot/massfunction-evolution.pdf")
p.clf()

rho0i=aa.critical_density0.to(uu.solMass/(uu.megaparsec)**3)
rho0=1.51*10**9 * 3840**3 / 10**9 # Modot/Mpc^3

js=[0,1,3,10,30,50]
for jj in js:
	print jj
	snList[jj]
	snL=glob.glob("/Volumes/data/BigMD/1Gpc_3840_Planck1/MFMC/*"+ str(snList[jj])+ "*MF.hist.dat")
	nnC=n.empty([6,799])
	nnS=n.empty([6,799])
	for iii in range(len(snL)):
		mmin,mmax,nnC[iii],nnS[iii]=n.loadtxt(snL[iii],unpack=True)

	nC=n.sum(nnC,axis=0)
	nS=n.sum(nnS,axis=0)
	mass=(mmin+mmax)/2.
	dLogM=n.median(mmax-mmin)
	sel=(nC>=10)&(nC<3840**3)
	p.plot(mass[sel], 10**(2*mass[sel])/rho0*nC[sel]/ dLogM/ volume,ls='None', marker='+', label= 'halos z='+str(numtoZ(snList[jj])))
	#p.plot(mass[sel],nS[sel]/dLogM/volume,'b+',label='sat')
	
p.grid()
p.xlabel('log halo mass')
p.ylabel(r'M$^2/\rho_0$ dN/dlogM ')
p.yscale('log')
#p.ylim((1e-7,1))
p.xlim((10,16))
p.legend(fontsize=9)
p.savefig("1Gpc_3840_Planck1/MFplot/massfunction-M2-evolution.pdf")
p.clf()


sys.exit()

import numpy as n
import cPickle
massB=n.arange(8,16,0.01)
vcirB=n.arange(0,4.5,0.01)
concB=n.arange(0,200,0.5)
for ii in range(len(snL)):
	print snL[ii]
	mtot, rvir, vcir, conc,distinct=n.loadtxt(snL[ii],unpack=True)
	cen=(distinct==0)
	sat=(cen==False)

	nnS,bb=n.histogram(n.log10(mtot[sat]),bins=massB)
	nnC,bb=n.histogram(n.log10(mtot[cen]),bins=massB)
	n.savetxt(snL[ii][:-6]+"MF.hist.dat",n.transpose([bb[:-1], bb[1:],nnC,nnS]))
	print "M"

	nnS,bb=n.histogram(n.log10(vcir[sat]),bins= vcirB)
	nnC,bb=n.histogram(n.log10(vcir[cen]),bins= vcirB)
	n.savetxt(snL[ii][:-6]+"VCIR.hist.dat",n.transpose([bb[:-1], bb[1:],nnC,nnS]))
	print"V"

	dataC=n.histogram2d(n.log10(mtot[cen]),conc[cen],bins=[massB,concB])
	dataS=n.histogram2d(n.log10(mtot[sat]),conc[sat],bins=[massB,concB])
	f=open(snL[ii][:-6]+"MCr.2dhist.pkl",'w')
	cPickle.dump([dataC,dataS],f)
	f.close()
	print "MC"

	dataC=n.histogram2d(n.log10(mtot[cen]),n.log10(vcir[cen]),bins=[massB, vcirB])
	dataS=n.histogram2d(n.log10(mtot[sat]),n.log10(vcir[sat]),bins=[massB, vcirB])
	f=open(snL[ii][:-6]+"MVr.2dhisti.pkl",'w')
	cPickle.dump([dataC,dataS],f)
	f.close()
	print "MV"

