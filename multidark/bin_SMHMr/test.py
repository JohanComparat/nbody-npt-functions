import XrayLuminosity
import numpy as n
import sys
import matplotlib.pyplot as p

xr = XrayLuminosity.XrayLuminosity()

print ( xr.psi(11., 1.1, n.arange(32,35,0.5)))
logM = 11
z = 1.1
pd = lambda x : n.piecewise(x, [x<32, x>=32], [lambda x : 0., lambda x : xr.psi(logM, z, x)] )
from scipy.integrate import quad
norm = quad(pd,32,35)[0]
pdn = lambda x : pd(x)/norm

x = n.arange(29,37, 0.1)
y = pdn(x)

p.plot(x,y)
p.yscale('log')
p.show()

from scipy import stats
class deterministic_gen(stats.rv_continuous):
	def _pdf(self, x):
		return pdn(x)


aa = deterministic_gen()
aa.rvs(size=1)

sys.exit()

import StellarMass
import numpy as n
from scipy.stats import norm
import glob
sm = StellarMass.StellarMass()
mhs = n.logspace(10,15,99)

ratio = sm.SMHMr(mhs,0.)
stellar_mass = sm.meanSM(mhs,0.)

print sm.sample_Ms( mhs, 0., scatter = 0.15 )
"""
import matplotlib.pyplot as p

for i in range(1000):
  #sms = 10**(sm.sample_Ms( mhs, 0., sigma_star = 0.15, loc = 1., scale = n.e**sm.meanSM(mhs, 0.) ))
  sms = 10**norm.rvs( loc = sm.meanSM(mhs, 0.), scale = 0.15 )
  #sms = 10**(sm.sample_Ms(mhs,0.15,sigma_star = 0.0225, scale=1.)-0.0225)
  p.loglog(mhs, sms/mhs, 'bx', alpha=0.2)

p.loglog(mhs, 10**stellar_mass/mhs, 'c--', lw=2, label='z=0')

mhs = n.logspace(10,15,98)
stellar_mass_1 = sm.meanSM(mhs,1.)

for i in range(1000):
  #sms = 10**(sm.sample_Ms( mhs, 0., sigma_star = 0.15, loc = 1., scale = n.e**sm.meanSM(mhs, 0.) ))
  sms = 10**norm.rvs( loc = sm.meanSM(mhs, 1.), scale = 0.15 )
  #sms = 10**(sm.sample_Ms(mhs,0.15,sigma_star = 0.0225, scale=1.)-0.0225)
  p.loglog(mhs, sms/mhs, 'r+', alpha=0.1)

p.loglog(mhs, 10**stellar_mass_1/mhs, 'm--', lw=2, label='z=1')

mhs = n.logspace(10,15,97)
stellar_mass_1 = sm.meanSM(mhs,2.)

for i in range(1000):
  #sms = 10**(sm.sample_Ms( mhs, 0., sigma_star = 0.15, loc = 1., scale = n.e**sm.meanSM(mhs, 0.) ))
  sms = 10**norm.rvs( loc = sm.meanSM(mhs, 2.), scale = 0.15 )
  #sms = 10**(sm.sample_Ms(mhs,0.15,sigma_star = 0.0225, scale=1.)-0.0225)
  p.loglog(mhs, sms/mhs, 'g+', alpha=0.05)

p.loglog(mhs, 10**stellar_mass_1/mhs, 'k--', lw=2, label='z=2')

p.grid()
p.xlabel('halo mass')
p.ylabel('stellar mass / halo mass')
p.legend()
p.savefig("/home/comparat/Desktop/smhmr.jpg")
p.clf()
"""
import astropy.io.fits as fits
import os
import time

aexp = 0.74230
z = 1./0.74230 -1.

fileList = n.array(glob.glob(os.path.join(os.environ['MD04'], "snapshots", "hlist_" + str(aexp) + "_SAM_Nb_*.fits" )))
fileList.sort()

# loop on the files
ii=0
for fileName in fileList:
  t0=time.time()
  outFile = os.path.join(os.environ['MD04'], "catalogs", os.path.basename(fileName)[:-5] + ".Ms.fits")
  print outFile
  hd = fits.open(fileName)
  Nhalo=len(hd[1].data['mvir'])

  col0 = fits.Column(name='Mgal_mvir_Mo13',format='D', array = norm.rvs( loc = sm.meanSM(10**hd[1].data['mvir'], z), scale = 0.15 ) )
  col0 = fits.Column(name='Mgal_m200c_Mo13',format='D', array = norm.rvs( loc = sm.meanSM(10**hd[1].data['m200c'], z), scale = 0.15 ) )
  col1 = fits.Column(name='Lx_agn',format='D', array = n.ones(Nhalo) )
  col2 = fits.Column(name='Lx_cluster',format='D', array = n.ones(Nhalo) )

  #define the table hdu 
  colArray = []
  for col in hd[1].columns :
      colArray.append(col)
  
  colArray.append(col0)
  colArray.append(col1)
  colArray.append(col2)

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
