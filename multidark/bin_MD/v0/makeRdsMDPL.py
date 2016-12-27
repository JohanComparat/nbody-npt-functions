import glob
import numpy as n
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
from scipy.interpolate import interp1d


rdList=["/data2/DATA/eBOSS/ELG/HAM-mocks/random-LC-4M.cat"]
catList=n.array(glob.glob("/data2/DATA/eBOSS/ELG/HAM-mocks-decam240/SHAM*_radecz.cat"))

print len(catList)

factor=10.
for ii in range(len(catList)):
	outname=catList[ii]+".random.cat"
	print outname
	ra,dec,z=n.loadtxt(catList[ii],unpack=True,usecols=(0,1,2))
	raR,decR=n.loadtxt(rdList[0],unpack=True,usecols=(0,1))

	dz=0.025
	z1=n.arange(0.5,1.11,dz)
	nn,bb,pp=p.hist(z,bins=z1)
	print nn.sum(),nn.sum(),nn.sum()*factor
	nFIN=int(nn.sum()*factor)+1
	p.clf()
	nz=interp1d((z1[1:]+z1[:-1])/2.,factor*nn)
	zs=n.arange(0.6,1.01,dz)
	rdsz=[]
	for i in range(len(zs)-1):
		inter=n.random.uniform(low=zs[i], high=zs[i+1], size=int(2* nz( zs[i]+dz/2. )))
		rdsz.append(inter)

	rds=n.hstack((rdsz))
	n.random.shuffle(rds)
	selRDS=(n.random.rand(len(raR))<float(nFIN)/len(raR))
	RR=rds[:len(raR[selRDS])]
	print "N final",len(raR[selRDS])
	n.savetxt(outname,n.transpose([raR[selRDS],decR[selRDS],RR]),fmt='%.8f %.8f %.4f')
	raR,decR,RR=0,0,0



