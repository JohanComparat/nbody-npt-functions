import sys
import glob
import os
import numpy as n
import astropy.io.fits as fits

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)

omega = lambda zz: cosmoMD.Om0*(1+zz)**3. / cosmoMD.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)

snList= n.array(glob.glob(os.path.join(os.environ["HOME"], 'MD04', "hlists", "hlist_*.list")))
snList.sort()
	
out_name = os.path.join(os.environ["HOME"], 'MD04', 'output_MD_0.4Gpc.fits')

names = n.array([el.split('_')[-1][:-5] for el in snList])
aexps =names.astype('float')

redshift = 1./n.array(aexps)-1.
dCom = cosmoMD.comoving_distance(redshift)
ids = n.argsort(redshift)
col0 = fits.Column(name='snap_name'	,format='A4', array = n.array(names)[ids] )
col2 = fits.Column(name='aexp'		,format='D', array = n.array(aexps)[ids] )
col3 = fits.Column(name='redshift'	,format='D', array = redshift[ids] )
col4 = fits.Column(name='comoving_distance',format='D', array= dCom.value[ids] )
age_yr = 10**9* cosmoMD.age(redshift).value
col5 = fits.Column(name='age_yr'		,format='D', array = age_yr[ids] )
array = cosmoMD.arcsec_per_kpc_comoving(redshift).value*3.6
col6 = fits.Column(name='deg_per_Mpc_comoving', format='D', array = array[ids] )


#define the table hdu 
hdu_cols  = fits.ColDefs([col0, col2, col3, col4, col5, col6])#, col7, col8, col9, col10])
tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
#define the header
prihdr = fits.Header()
prihdr['sim'] = 'SMD'
prihdr['author'] = 'JC'
prihdu = fits.PrimaryHDU(header=prihdr)
#writes the file
thdulist = fits.HDUList([prihdu, tb_hdu])
#os.system("rm "+out_name)
thdulist.writeto(out_name)

rho_crit = cosmoMD.critical_density(redshift[ids]).to(u.Msun/u.Mpc**3.)
delta_vir = DeltaVir_bn98(redshift[ids])

n.savetxt("/u/joco/data/MD/MD_0.4Gpc/hlists_MD_0.4Gpc.ascii", n.transpose([aexps[ids], redshift[ids], age_yr[ids], rho_crit.value, delta_vir]))
