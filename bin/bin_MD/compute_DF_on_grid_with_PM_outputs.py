#  cd pySU/pyMultidark/trunk/bin

import numpy as n
import os
from os.path import join
from astropy.io import fits
import time
import fortranfile
import cPickle
from scipy.io import FortranFile as ff

DFdir = "/home2/jcomparat/eBOSS-LC/PM_runs/PM_C/PM_L650_g1612_Dgt120_am100_Dlt240/"

inFiles = n.array(["PMcrd.DAT","PMcrs1.0136.DAT"])

bins = n.hstack((-1,n.arange(0,10,0.1),n.arange(10,100,10),n.arange(100,1000,100),n.arange(1000,10000,1000)))


#for infi in inFiles:

infi = inFiles[1]
DFfile = join(DFdir,infi)
f=ff(DFfile, 'r')
record = f.read_record([('x', 'f4'),('y', 'f4'),('z', 'f4'),('vx', 'f4'),('vy', 'f4'),('vz', 'f4')])# ('b', '<i4', (3,3))])

f = open(DFfile, 'rb')
f.close()

f = fortranfile.FortranFile(DFfile, endian='=', header_prec='i')

2i7,2a,3x,i9
10x,a,i5,a,4i11

i5 : integer 5 positions
3x : 3 horizontal positionning space ?
2a : 2 characters


out = f.readRecord()
xlines = f.readlines()


gridx, gridy, gridz = f.readInts()
res = n.empty((gridx, len(bins)-1))
for kk in range(gridx):
	DF = f.readReals()
	res[kk],bi = n.histogram(DF,bins=bins)

f.close()
result = n.sum(res,axis=0)
path_to_outputCat =  join(mockDir,infi[:-4] + "_DFhist_linearBins.dat")
f=open(path_to_outputCat, 'w')
cPickle.dump([bins,result],f)
f.close()
