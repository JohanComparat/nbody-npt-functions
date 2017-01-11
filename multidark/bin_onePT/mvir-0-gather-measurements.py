import glob
from os.path import join
import numpy as n
import astropy.io.fits as fits
import lib_functions_1pt as lib
import os
import sys

#Quantity studied
version = 'v3'
qty = "mvir"


# one point function lists
fileC = n.array(glob.glob( join(os.environ['MD_DIR'], "MD_*Gpc*",  version, qty,"out_*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['MD_DIR'], "MD_*Gpc*", version, qty,"out_*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['MD_DIR'], "MD_*Gpc*", version, qty,"out_*_Satellite_JKresampling.pkl")))

fileC.sort()
fileS.sort()
fileB.sort()
print len(fileC), len(fileB), len(fileS)

print "considers ",len(fileC), qty , " function files"

for ii, el in enumerate(fileC):
	print el
	print fileS[ii]
	print fileB[ii]
	lib.convert_pkl_mass(fileC[ii], fileS[ii], fileB[ii], qty)


fileC = n.array(glob.glob( join(os.environ['DS_DIR'], version, qty,"ds*_Central_JKresampling.pkl")))
fileB = n.array(glob.glob( join( os.environ['DS_DIR'], version, qty,"ds*_"+qty+"_JKresampling.bins")))
fileS = n.array(glob.glob( join( os.environ['DS_DIR'], version, qty,"ds*_Satellite_JKresampling.pkl")))

fileC.sort()
fileS.sort()
fileB.sort()
print len(fileC), len(fileB), len(fileS)

print "considers ",len(fileC), qty , " function files"

for ii, el in enumerate(fileC):
	print el
	print fileS[ii]
	print fileB[ii]
	lib.convert_pkl_mass(fileC[ii], fileS[ii], fileB[ii], qty)
	

print qty
af = n.array(glob.glob(join(os.environ['MVIR_DIR'], "data", "*_"+qty+".fits") ) )
print af[0]
d0 = fits.open(af[0])[1].data
#print len(d0['log_mvir']), d0['log_mvir']

for ii in range(1,len(af),1):
	d1 = fits.open(af[ii])[1].data
	d0 = n.hstack((d0,d1))

hdu2 = fits.BinTableHDU.from_columns(d0)

writeName = join(os.environ['MVIR_DIR'], qty+"_summary.fits")
if os.path.isfile(writeName):
	os.remove(writeName)
	
hdu2.writeto( writeName )


"""
sys.exit()
# rebinning here
#solve bins = 0 problem

n.arange()
n.hstack((n.arange(8,14,0.25), n.arange(14,16,0.05)))

#if logmvir < 14 :
Nrb = 5.
idSplit = int(n.searchsorted(d0['log_mvir'],14)/Nrb)*Nrb
split_array= lambda array: [array[:idSplit], array[idSplit:]]

#variables :
# 
def rebinMe(trb, mod, Nrb = 5):
	# split
	part1, part2 = split_array(trb)
	# rebin
	take_middle_val = lambda part: part[2::Nrb]
	take_mean_val = lambda part: (part[0::Nrb] + part[1::Nrb] + part[2::Nrb] + part[3::Nrb] + part[4::Nrb])/Nrb.
	take_sum_val = lambda part: part[0::Nrb] + part[1::Nrb] + part[2::3] + part[3::Nrb] + part[4::Nrb]
	if mode == 'middle' :
		part1b = take_middle_val(part1)
	if mode == 'mean' :
		part1b = take_mean_val(part1)
	if mode == 'sum' :
		part1b = take_sum_val(part1)
	return n.hstack((part1b, part2))

trb = d0['log_mvir']
mode = 'middle'
trb_o = rebinMe(trb, mode)
"""

