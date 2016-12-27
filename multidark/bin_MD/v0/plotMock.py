import random
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import astropy.io.fits as pf
import numpy as n
from scipy.interpolate import interp1d
import scipy.stats as st
import glob
mockDir="/data2/DATA/eBOSS/ELG/HAM-mocks-decam240/"

print "outputs stored : ",mockDir

mockList=glob.glob(mockDir+"*_allCols.cat.gz")
for mockFile in mockList:
	ra,dec,z,mass,cs=n.loadtxt(mockFile,unpack=True)
	cen=(cs==-1)
	sat=(cen==False)
	f_cen=len(ra[cen])/float(len(ra))
	f_sat=1-f_cen
	p.figure(0,(5,5))
	p.axes([0.17,0.17,0.78,0.78])
	p.plot(z,n.log10(mass),'r+',alpha=0.1,rasterized=True,label="f sat "+str(n.round(100*f_sat,2))+"%")
	#p.plot(z[sat],n.log10(mass[sat]),'b+',rasterized=True,label="f sat "+str(n.round(100*f_sat,2))+"%")
	p.legend()
	p.xlabel('redshift')
	p.ylabel(r'$\log(M/M_\odot)$')
	p.savefig(mockFile[:-7]+".pdf")
	p.clf()

