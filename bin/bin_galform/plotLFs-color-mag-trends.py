#! /usr/bin/env python

"""
This script produces quality plots to check that the LFs are fine compared to simumlations.
"""

import sys
import os
from os.path import join
data_dir = os.environ['DATA_DIR']
import glob

from lib_plot import *
#from lineListAir import *
SNlim = 5

# "D:\data\LF-O\LFmodels\data\trends_color_mag\O2_3728-VVDSDEEPI24-z0.947.txt"

plotDir="/home/comparat/database/Simulations/galform-lightcone/products/emissionLineLuminosityFunctions/plots/"

dir="/home/comparat/database/Simulations/galform-lightcone/products/emissionLineLuminosityFunctions/O2_3728/"

lf_measurement_files_ref=n.array(glob.glob(dir+"*MagLimR-24.2-z0.7*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*MagLimR-*z0.7*.txt"))
lf_measurement_files.sort()


dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)
phiRatio = n.empty([ len(lf_measurement_files), len(dataRef[0]) ])
label = n.array([ "R<23.0","R<23.5", "R<24.2"])
for ii,el in  enumerate(lf_measurement_files) :
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	
	
imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((1e40,1e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=2)
p.savefig(join(plotDir,"trends_O2_3728_Rmag-z0.7.pdf"))
p.clf()


lf_measurement_files_ref=n.array(glob.glob(dir+"*MagLimR-24.2-z0.9*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*MagLimR-*z0.9*.txt"))
lf_measurement_files.sort()


dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)
phiRatio = n.empty([ len(lf_measurement_files), len(dataRef[0]) ])
label = n.array([ "R<23.0","R<23.5", "R<24.2"])
for ii,el in  enumerate(lf_measurement_files) :
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	
	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((1e40,1e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=2)
p.savefig(join(plotDir,"trends_O2_3728_Rmag-z0.9.pdf"))
p.clf()



lf_measurement_files_ref=n.array(glob.glob(dir+"*MagLimR-24.2-z1.*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*MagLimR-*z1.*.txt"))
lf_measurement_files.sort()


dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)
phiRatio = n.empty([ len(lf_measurement_files), len(dataRef[0]) ])
label = n.array([ "R<23.0","R<23.5", "R<24.2"])
for ii,el in  enumerate(lf_measurement_files) :
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	
	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((1e40,1e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=2)
p.savefig(join(plotDir,"trends_O2_3728_Rmag-z1.2.pdf"))
p.clf()






########################################33


lf_measurement_files_ref=n.array(glob.glob(dir+"*MagLimI-24-z1.*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*MagLimI-*z1.*.txt"))
lf_measurement_files.sort()


dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)
phiRatio = n.empty([ len(lf_measurement_files), len(dataRef[0]) ])
label = n.array([ "I<22.5", "I<23.0","I<23.5","I<24.0"])
for jj,el in  enumerate(lf_measurement_files) :
	data= n.loadtxt( el, unpack=True)
	phiRatio[jj] = data[3] / dataRef[3]	
	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((1e40,1e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=2)
p.savefig(join(plotDir,"trends_O2_3728_Imag-z1.2.pdf"))
p.clf()



lf_measurement_files_ref=n.array(glob.glob(dir+"*MagLimI-24-z0.9*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*MagLimI-*z0.9*.txt"))
lf_measurement_files.sort()


dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)
phiRatio = n.empty([ len(lf_measurement_files), len(dataRef[0]) ])
label = n.array([ "I<22.5", "I<23.0","I<23.5","I<24.0"])
for jj,el in  enumerate(lf_measurement_files) :
	data= n.loadtxt( el, unpack=True)
	phiRatio[jj] = data[3] / dataRef[3]	
	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((1e40,1e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=2)
p.savefig(join(plotDir,"trends_O2_3728_Imag-z0.9.pdf"))
p.clf()



lf_measurement_files_ref=n.array(glob.glob(dir+"*MagLimI-24-z0.7*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*MagLimI-*z0.7*.txt"))
lf_measurement_files.sort()


dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)
phiRatio = n.empty([ len(lf_measurement_files), len(dataRef[0]) ])
label = n.array([ "I<22.5", "I<23.0","I<23.5","I<24.0"])
for jj,el in  enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[jj] = data[3] / dataRef[3]	
	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((1e40,1e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=2)
p.savefig(join(plotDir,"trends_O2_3728_Imag-z0.75.pdf"))
p.clf()



#####################################3
#####################################3
# R-Z
#####################################3
#####################################3

lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSrz_gt_0.0-z0.7*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSrz_?t_*z0.7*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

label = n.array(["r-z>0", "r-z>0.5", "r-z>1", "r-z>1.5", "r-z<1", "r-z<1.5", "r-z<2"])
phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_RZ-z0.75.pdf"))
p.clf()


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSrz_gt_0.0-z0.9*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSrz_?t_*z0.9*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_RZ-z0.9.pdf"))
p.clf()


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSrz_gt_0.0-z1.*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSrz_?t_*z1.*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_RZ-z1.2.pdf"))
p.clf()



#####################################3
#####################################3
# G-R
#####################################3
#####################################3


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSgr_gt_0.0-z0.7*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSgr_?t_*z0.7*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

label = n.array(["g-r>0", "g-r>0.5", "g-r>1", "g-r>1.5", "g-r<1", "g-r<1.5", "g-r<2"])
phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_GR-z0.75.pdf"))
p.clf()


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSgr_gt_0.0-z0.9*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSgr_?t_*z0.9*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

label = n.array(["g-r>0", "g-r>0.5", "g-r>1", "g-r>1.5", "g-r<1", "g-r<1.5", "g-r<2"])
phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_GR-z0.9.pdf"))
p.clf()


lf_measurement_files_ref=n.array(glob.glob(dir+"*VVDSgr_gt_0.0-z1.*.txt"))
lf_measurement_files=n.array(glob.glob(dir+"*VVDSgr_?t_*z1.*.txt"))
lf_measurement_files.sort()

dataRef = n.loadtxt( lf_measurement_files_ref[0], unpack=True)

label = n.array(["g-r>0", "g-r>0.5", "g-r>1", "g-r>1.5", "g-r<1", "g-r<1.5", "g-r<2"])
phiRatio = n.empty([ 7, len(dataRef[0]) ])
for ii, el in enumerate(lf_measurement_files):
	data= n.loadtxt( el, unpack=True)
	phiRatio[ii] = data[3] / dataRef[3]	

imin = n.argmax(dataRef[6])-1
p.figure(0,(6,6))
for jj in range(len(label)):
	p.plot(dataRef[2][imin:],phiRatio[jj][imin:],label=label[jj])

p.xlabel(r'$log_{10}(L[O_{II}])$ [erg s$^{-1}$]')
p.ylabel(r'$\Phi/\Phi_{ref}$')
p.xscale('log')
p.xlim((7e40,5e43))
p.ylim((-0.05,1.05))
p.grid()
p.legend(loc=4)
p.savefig(join(plotDir,"trends_O2_3728_I22.5_GR-z1.2.pdf"))
p.clf()

