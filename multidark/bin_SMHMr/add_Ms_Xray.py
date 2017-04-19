import StellarMass
import XrayLuminosity

import numpy as n
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d

import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

print " set up box, and redshift "
aexp = 0.74230
z = 1./0.74230 -1.
fileList = n.array(glob.glob(os.path.join(os.environ['MD04'], "snapshots", "hlist_" + str(aexp) + "_SAM_Nb_*.fits" )))
fileList.sort()
print fileList
# set up the stellar mass computation
sm = StellarMass.StellarMass()
mhs = n.logspace(10,15,99)

ratio = sm.SMHMr(mhs,0.)
stellar_mass = sm.meanSM(mhs,0.)

# set up the x ray lambda SAR
xr = XrayLuminosity.XrayLuminosity()
logMs = n.arange(6.5,12.5,0.01)
cdfs_interpolations = []
XXS = n.arange(32,36.1,0.1)
for jj,mass in enumerate(logMs):
	pd = lambda ll : xr.psi(ll, logM=mass, z=z)
	norm = quad( pd, 32, 36)[0]
	cdfs_interpolations.append( interp1d(n.array([quad( pd, 32, X)[0] for X in XXS ])/norm,XXS) )

cdfs_interpolations = n.array(cdfs_interpolations)

print " loop on the files "
ii=0
for fileName in fileList:
  t0=time.time()
  outFile = os.path.join(os.environ['MD04'], "catalogs", os.path.basename(fileName)[:-5] + ".Ms.fits")
  print outFile
  hd = fits.open(fileName)
  Nhalo=len(hd[1].data['mvir'])
  Mgal_mvir_Mo13 = norm.rvs( loc = sm.meanSM(10**hd[1].data['mvir'], z), scale = 0.15 )
  randomX = n.random.rand(len(Mgal_mvir_Mo13))
  indexes = n.searchsorted(logMs,Mgal_mvir_Mo13)
  lambda_sar_Bo16 = n.array([ cdfs_interpolations[indexes[ii]](randomX[ii]) for ii in range(Nhalo) ])
  col0 = fits.Column(name='Mgal_mvir_Mo13',format='D', array = Mgal_mvir_Mo13 )
  col1 = fits.Column(name='Mgal_m200c_Mo13',format='D', array = norm.rvs( loc = sm.meanSM(10**hd[1].data['m200c'], z), scale = 0.15 ) )
  col2 = fits.Column(name='lambda_sar_Bo16',format='D', array = lambda_sar_Bo16 )
  col3 = fits.Column(name='Lx_cluster',format='D', array = n.ones(Nhalo) )

  #define the table hdu 
  colArray = []
  for col in hd[1].columns :
      colArray.append(col)
  
  colArray.append(col0)
  colArray.append(col1)
  colArray.append(col2)
  colArray.append(col3)
  
  hdu_cols  = fits.ColDefs(colArray)
  tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
  #define the header
  prihdr = fits.Header()
  prihdr['author'] = 'JC'
  prihdu = fits.PrimaryHDU(header=prihdr)
  #writes the file
  thdulist = fits.HDUList([prihdu, tb_hdu])
  #os.system("rm "+self.snl[ii][:-5]+"_PM_Nb_"+str(Nb)+".fits")
  if os.path.isfile(outFile):
    os.system("rm "+outFile)
    
  thdulist.writeto(outFile)
  print time.time()-t0


