import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

# specific functions
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d

summ = fits.open(os.path.join(os.environ["MD10"], 'output_MD_1.0Gpc.fits'))[1].data	
volume = (1000/0.6777)**3.

zmin, zmax, zmean, dN_dVdz_lrg, dN_dVdz_elg, dN_dVdz_qso = n.loadtxt(os.path.join(os.environ['GIT_NBODY_NPT'], "data/NZ/nz_8.txt"), unpack=True)

dN_dz_elg = interp1d(n.hstack((-1., zmin[0], zmean, zmax[-1], 3.)), volume * n.hstack((0., 0., dN_dVdz_elg, 0., 0.)) )
dN_dz_lrg = interp1d(n.hstack((-1., zmin[0], zmean, zmax[-1], 3.)), volume * n.hstack((0., 0., dN_dVdz_lrg, 0., 0.)) )
dN_dz_qso = interp1d(n.hstack((-1., zmin[0], zmean, zmax[-1], 3.)), volume * n.hstack((0., 0., dN_dVdz_qso, 0., 0.)) )

# determine dz
z_middle = (summ['redshift'][1:]+summ['redshift'][:-1])*0.5
z_mins = n.hstack((summ['redshift'][0], z_middle))
z_maxs = n.hstack((z_middle, summ['redshift'][-1]))
z_snap = summ['redshift']
dz = z_maxs - z_mins

N_elg_per_snap = (dN_dz_elg(z_snap) * dz).astype('int')
N_lrg_per_snap = (dN_dz_lrg(z_snap) * dz).astype('int')
N_qso_per_snap = (dN_dz_qso(z_snap) * dz).astype('int')

def SHAM_4MOST(fileList_snap, z, dz, nTargets =[1., 1., 1.]):

for ii, el in enumerate(summ):
	nTargets = dz[ii]*n.array([dN_dz_lrg(z), dN_dz_elg(z), dN_dz_qso(z)])
	fileList_snap = n.array(glob.glob(os.path.join(os.environ["MD10"], 'work_agn', 'out_'+el['snap_name']+'_SAM_Nb_?.fits')))
	fileList_snap.sort()
	print( el )
	print( fileList_snap )
	print( nTargets )
	SHAM_4MOST(fileList_snap, z, dz[ii], nTargets =[1., 1., 1.])


from populateHaloCatalogsLib import *
import os


interp1d()

p_elg=[10**(12.2),0.25]
p_qso=[10**(12.7),0.25]
p_lrg=[10**(12.2),0.25]

ids=[[],[],[], [],[],[]]
oN=["","","", "","",""]
for ii in range(len(zmin)):
	nG=int(nGal_Deg2[ii]* area)+1
	print "gets all halos for ",zmin[ii],"<z<",zmax[ii], "with col5 to mock ", nG, " galaxies." 
	IDhz_c,QTY_c,nn_c,bb_c=get_distrib_QTY_cen(hdu, colN='col5', zmin=zmin[ii], zmax=zmax[ii])
	IDhz_s,QTY_s,nn_s,bb_s=get_distrib_QTY_sat(hdu, colN='col5', zmin=zmin[ii], zmax=zmax[ii])
	for jj in range(len(ps)):
		p1,p2,p3=ps[jj] #4e13,1e13,0.5111
		oN[jj]="SHAM_norm-mean"+str(p1)+"-sig"+str(p2)+"-fsat"+str(p3)+"_ELGdensity"
		ids[jj].append(selectGaussian_fsat(p1,p2,p3, nG,IDhz_c, QTY_c,IDhz_s, QTY_s ))

for jj in range(len(ps)):
	writerCatsAll(oN[jj],n.sort(n.hstack((ids[jj]))))
