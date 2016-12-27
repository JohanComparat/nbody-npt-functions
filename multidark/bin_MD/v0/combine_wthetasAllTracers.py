import numpy as n
import glob

haloList=glob.glob("/data2/DATA/eBOSS/all-tracer-mocks/SHAM_*_radecz.cat")
for el in haloList:
	dataPath=el
	xx0,yy0,y0E=n.loadtxt(dataPath[:-4]+"_angular_d3.dat",unpack=True,usecols=(0,1,2))
	xx1,yy1,y1E=n.loadtxt(dataPath[:-4]+"_angular_d2.dat",unpack=True,usecols=(0,1,2))
	xx2,yy2,y2E=n.loadtxt(dataPath[:-4]+"_angular_d1.dat",unpack=True,usecols=(0,1,2))
	xx=n.hstack((xx0,xx1[1:],xx2[1:]))
	yy=n.hstack((yy0,yy1[1:],yy2[1:]))
	yE=n.hstack((y0E,y1E[1:],y2E[1:]))
	xE=n.hstack((n.ones_like(xx0)*xx0[0],n.ones_like(xx1[1:])*xx1[0], n.ones_like(xx2[1:])*xx2[0]))
	n.savetxt(dataPath[:-4]+"_"+"angular_combined.dat", n.transpose([xx,yy,yE,xE]))

