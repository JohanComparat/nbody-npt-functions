from HaloSelection import *


"""
Library of function to create halo catalogs matched to a density.

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

"""
import random
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import astropy.io.fits as fits
import numpy as n
from scipy.interpolate import interp1d
import scipy.stats as st
import os
from os.path import join

# initialize data side

# what tracer we are dealing with :
tracer_dir = join(os.environ['DATA_DIR'],"ELG")

# where is the NZ :
nz_dir = join(tracer_dir, "observations/NZ")
NZfile = join(nz_dir, "nz-fisherGRIW1.dat")

# loads the NZ, needs to be per deg2
zminIN, zmaxIN, nGal_Deg2IN = n.loadtxt( NZfile, unpack = True, usecols = (0,1,2) )
ok = (nGal_Deg2IN>0) & (zmaxIN<1.25) & (zminIN > 0.39)
zmin, zmax, nGal_Deg2  =  zminIN[ok], zmaxIN[ok], nGal_Deg2IN[ok]

# where the outputs will be stored :
mockOutpuName = "mocks_fischerGRIW1"
mockOutput_dir = join(tracer_dir,mockOutpuName)
os.system('mkdir ' + mockOutput_dir)
mockName = "tryMocks-sham"
print "outputs stored : ",mockOutput_dir
print "N(z) from ", NZfile


# initialize the lightcone to extract the mock from
lcDir = "/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/lightcones/lc_square_0.1z1.4/"
lcFile = join(lcDir,"lightcone_MD_2.5Gpc_0.4z1.4.fits")

# loads the LC
hdu = fits.open(lcFile)
print lcFile, " loaded, columns:"
print hdu[1].data.dtype
# [('ra', '>f8'), ('dec', '>f8'), ('z_real_space', '>f8'), ('z_redshift_space', '>f8'), ('v_parallel', '>f8'), ('id', '>i8'), ('num_prog', '>i2'), ('pid', '>i8'), ('upid', '>i8'), ('mvir', '>f4'), ('rvir', '>f8'), ('rs', '>f8'), ('vrms', '>f8'), ('vmax', '>f8'), ('Spin', '>f4'), ('M200b', '>f8'), ('M200c', '>f8'), ('b_to_a', '>f8'), ('c_to_a', '>f8'), ('Halfmass_Radius', '>f4'), ('Macc', '>f4'), ('Mpeak', '>f4'), ('Vacc', '>f8'), ('Vpeak', '>f8'), ('Acc_Rate_Inst', '>f4'), ('Acc_Rate_100Myr', '>f4'), ('Acc_Rate_Mpeak', '>f4')]

# properties of the lightcone
area = (2*30.)**2.

# boolean arrays that discriminate central and sat halos
cen = (hdu[1].data['pid'] ==  -1)
sat = (cen ==  False)


MMD = MultiDarkMock(hdu, area, mockOutput_dir, mockName, zmin, zmax, nGal_Deg2)

MMD.mockName = "tryMocks-sham"
MMD.make_sham_catalog(colN='mvir')
MMD.write_full_catalog_fits()
MMD.write_catalog_ascii()

MMD.mockName = "tryMocks-shamInc"
MMD.make_shamIncomplete_catalog( colN = 'mvir', incompletenessFactor = 0.5 * n.ones_like(MMD.zmin))
MMD.write_full_catalog_fits()

MMD.mockName = "tryMocks-shamMax"
MMD.make_shamMAX_catalog( colN='mvir', maxQTY = n.ones_like(MMD.zmin) * 1e13 )
MMD.write_full_catalog_fits()

p1 = 12.8
p2 = 12.4
MMD.mockName = "tryMocks-gaussian_mean_"+str(p1)+"_sig_"+str(p2)
make_catalog_Gaussian( colN='mvir', means = n.ones_like(MMD.zmin) * 10**p1, scatters = n.ones_like(MMD.zmin) * 10**p2 )
MMD.write_full_catalog_fits()

p1 = 12.8
p2 = 12.4
p3 = 0.20
MMD.mockName = "tryMocks-gaussian_mean_"+str(p1)+"_sig_"+str(p2)+"_fsat_"+str(p3)
MMD.make_catalog_GaussianFsat( colN='mvir', means = n.ones_like(MMD.zmin) * 10**p1, scatters = n.ones_like(MMD.zmin) * 10**p2 , fsats = n.ones_like(MMD.zmin) * p3 )
MMD.write_full_catalog_fits()

p1 = 12.8
p2 = 12.4
MMD.mockName = "tryMocks-lognorm_mean_"+str(p1)+"_sig_"+str(p2)
make_LogNorm_catalog( colN='mvir', means = n.ones_like(MMD.zmin) * 10**p1, scatters = n.ones_like(MMD.zmin) * 10**p2 )
MMD.write_full_catalog_fits()


# load lightcone

# loop over mocks

sats=n.arange(0.,0.31,0.025)
mMean=n.array([4,8,16,32,64,128])*10**(12)
ps=[]
ids,oN=[],[]
for el in mMean:
	for fs in sats:
		ps.append([el,el,fs])
		ids.append([])
		oN.append("")
		ps.append([el,el/2.,fs])
		ids.append([])
		oN.append("")
		ps.append([el,el/4.,fs])
		ids.append([])
		oN.append("")

print len(ps), len(ids), len(oN)

for ii in range(len(zmin)):
	nG=int(nGal_Deg2[ii]* area)+1
	print "gets all halos for ",zmin[ii],"<z<",zmax[ii], "with col5 to mock ", nG, " galaxies." 
	IDhz_c,QTY_c,nn_c,bb_c=get_distrib_QTY_cen(hdu, colN='mvir', zmin=zmin[ii], zmax=zmax[ii])
	IDhz_s,QTY_s,nn_s,bb_s=get_distrib_QTY_sat(hdu, colN='mvir', zmin=zmin[ii], zmax=zmax[ii])
	for jj in range(len(ps)):
		p1,p2,p3=ps[jj] #4e13,1e13,0.5111
		oN[jj]="SHAM_norm-mean"+str(p1)+"-sig"+str(p2)+"-fsat"+str(p3)+"_ELGdensity"
		ids[jj].append(selectGaussian_fsat(p1,p2,p3, nG,IDhz_c, QTY_c,IDhz_s, QTY_s ))

for jj in range(len(ps)):
	writerCats(oN[jj],n.sort(n.hstack((ids[jj]))))


