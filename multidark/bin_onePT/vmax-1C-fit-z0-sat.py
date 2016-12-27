from os.path import join
import numpy as n
import astropy.io.fits as fits
import os
import lib_functions_1pt as lib
import cPickle
import astropy.cosmology as co
cosmo = co.Planck13
import glob
import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p

from scipy.interpolate import interp1d
from scipy.misc import derivative
#Quantity studied
#=================
qty = "vmax"

cos = 'sat'

# working directory
#=================
dir = join(os.environ['MULTIDARK_LIGHTCONE_DIR'], qty)
# loads summary file
data = fits.open( join(dir, "MD_"+qty+"_summary.fits"))[1].data
ZMIN = 0. 
ZMAX = 2.3
MD04=(data["boxName"]=='MD_0.4Gpc')
MD10=(data["boxName"]=='MD_1Gpc')
MD25=(data["boxName"]=='MD_2.5Gpc')
MD40=(data["boxName"]=='MD_4Gpc')
MD25NW=(data["boxName"]=='MD_2.5GpcNW')
MD40NW=(data["boxName"]=='MD_4GpcNW')
MD40z = (MD40) & (data['redshift']<ZMAX)& (data['redshift']>ZMIN)
zzs = set(data['redshift'][MD40z])

param_z0_file=open(join(dir, "vmax-"+cos+"-diff-function-z0-params.pkl"), 'r')
outCF = cPickle.load(param_z0_file)
param_z0_file.close()
A0, v0, a0, b0 = outCF[0]

# fitting function parameters
#=================
def fitZ(zmin, zmax, suffix):
	NminCount = 100
	# 1 percent deviation from mass function given the sigma 8 deviation
	limits_04 =  [125, 450] #max : 2e13#[100, 1000]
	limits_10 =  [250, 800] #max : 2e14
	limits_25 =  [600, 1100] #max : 5e14
	limits_40 =  [800, 1400] #max : 2e14

	p0 = n.array([A0, v0, a0, b0])

	#zmin = -0.01
	#zmax = 0.001

	#=================
	# DATA
	#=================

	# redshift selection
	zSel = lib.zSelection( data, zmin, zmax )
	# mass selection
	mSel = lib.vSelection(data, qty, limits_04, limits_10, limits_25,limits_40) 
	# minimum number counts selection
	nSel = lib.nSelection(data, NminCount, cos)
	# altogether
	ok = (zSel) & (mSel) & (nSel)
	# selection per box :


	# x coordinates definition
	#=================
	vmax = data[qty]
	log_vmax = n.log10(vmax)
	#print len(vmax), n.min(vmax), n.max(vmax)

	# y coordinates
	#=================
	norm = (100)**3. /(cosmo.H(data["redshift"]).value)**6.
	log_VF = n.log10( norm * vmax**3. * data["dNdVdlnV_"+cos])
	log_VF_c = n.log10( norm * vmax**3. * data["dNdVdlnV_"+cos+"_c"])
	#print n.min(log_VF), n.max(log_VF)

	# error on y position
	#=================
	error = (data["std90_pc_"+cos]**2. + data["dN_counts_"+cos]**(-1.))**(0.5)
	print "n points", len(log_vmax[ok])
	if len(log_vmax[ok])>50 :
		pOpt, pCov = lib.fit_vmax_function_z0(data[ok], x_data = log_vmax[ok], y_data = log_VF[ok], y_err = error[ok], p0 = p0, cos = cos, mode = "curve_fit", suffix = suffix)

zba = n.array(list(zzs))
zba.sort()
#nna = (zba[1:]+zba[:-1])/2.
#zbb = n.hstack(( zba[0]-0.1, nna, zba[-1]+0.1 ))

for ii, zz in enumerate(zba) :
	print '====================='
	print 'redshift bin', str(zz)
	#fitZ(zmin=zbb[ii], zmax=zbb[ii+1], suffix=str(zz))
	fitZ(zmin=zz-0.001, zmax=zz+0.001, suffix=str(zz))

print '====================='
print '====================='
print '====================='
print 'Z=0'
print '====================='
print '====================='
fitZ(zmin=-0.001, zmax=+0.001, suffix=str(0.0))
outfiles = n.array(glob.glob(join(dir,"vmax-"+cos+"-diff-function-params-*.pkl")))
outfiles.sort()

ps = []
psErr = []
ZZ = []
for outfile in outfiles:
	f=open(outfile, 'r')
	ZZ.append(float(os.path.basename(outfile).split('-')[-1][:-4]))
	pOpt, pCov = cPickle.load(f)
	ps.append(pOpt)
	psErr.append(pCov.diagonal()**0.5)

pA, pv0, palpha, pbeta = n.transpose(ps)
pA_e, pv0_e, palpha_e, pbeta_e = n.transpose(psErr)
ZZ = n.array(ZZ)

ZZr = n.arange(0., ZMAX+0.02, 0.01)

def plotPoly(ndeg):
	fig = p.figure(10, (8, 8))
	ax1 = fig.add_subplot(221)
	ax1.errorbar(ZZ, pA, yerr = pA_e, rasterized=True, fmt='none')
	ppp = n.polyfit(ZZ, pA, deg=ndeg, w= 1/pA_e)
	print 'A', ppp
	ax1.plot(ZZr, n.polyval(ppp, ZZr), label=ppp)
	ax1.set_xlabel('z')
	ax1.set_ylabel('A')
	ax1.set_ylim((-2,0.5))
	ax1.set_xlim((-0.1,ZMAX+0.1))
	ax1.grid()
	ax1.set_title(n.round(ppp,3), fontsize=9)
	ax2 = fig.add_subplot(222)
	ax2.errorbar(ZZ, pv0, yerr = pv0_e, rasterized=True, fmt='none')
	ppp = n.polyfit(ZZ, pv0, deg=ndeg, w= 1/pv0_e)
	print 'v0', ppp
	ax2.plot(ZZr, n.polyval(ppp, ZZr), label=ppp)
	ax2.set_xlabel('z')
	ax2.set_ylabel('v0')
	ax2.set_ylim((2,4))
	ax2.set_xlim((-0.1,ZMAX+0.1))
	ax2.grid()
	ax2.set_title(n.round(ppp,3), fontsize=9)
	ax3 = fig.add_subplot(223)
	ax3.errorbar(ZZ, palpha, yerr = palpha_e, rasterized=True, fmt='none')
	ppp = n.polyfit(ZZ, palpha, deg=ndeg, w= 1/palpha_e)
	print 'alpha', ppp
	ax3.plot(ZZr, n.polyval(ppp, ZZr), label=ppp)
	ax3.set_xlabel('z')
	ax3.set_ylabel('alpha')
	ax3.set_ylim((0.5,2.5))
	ax3.set_xlim((-0.1, ZMAX+0.1))
	ax3.set_title(n.round(ppp,3), fontsize=9)
	ax3.grid()
	ax4 = fig.add_subplot(224)
	ax4.errorbar(ZZ, pbeta, yerr = pbeta_e, rasterized=True, fmt='none')
	ppp = n.polyfit(ZZ, pbeta, deg=ndeg, w= 1/pbeta_e)
	print 'beta', ppp
	ax4.plot(ZZr, n.polyval(ppp, ZZr), label=ppp)
	ax4.set_title(n.round(ppp,3), fontsize=9)
	ax4.set_xlabel('z')
	ax4.set_ylabel('beta')
	ax4.set_ylim((-3,0))
	ax4.set_xlim((-0.1,ZMAX+0.1))
	ax4.grid()
	p.savefig(join(dir, "parameters-ztrend-"+cos+str(ndeg)+".png"))
	p.clf()


plotPoly(1)
plotPoly(2)
plotPoly(3)
plotPoly(4)
plotPoly(5)





