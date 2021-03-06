import MultiDark as MD
import sys
import glob
import os
import numpy as n
import astropy.io.fits as fits
import astropy.units as uu

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
omega = lambda zz: cosmoMD.Om0*(1+zz)**3. / cosmoMD.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)
	
snList= n.array(glob.glob(os.path.join(os.environ["MD10"], "hlists", "hlist_*.list")))
snList.sort()

def get_first(path_2_snapshot):
	fl = MD.fileinput.input(path_2_snapshot)
	aexp = float(os.path.basename(path_2_snapshot)[:-5].split('_')[1])
	l0 = fl.readline()[1:]
	fl.close()
	return l0.split(), aexp

names, colN, aexps = [], [], []
for snLi in snList:
	lin1, aaa = get_first(snLi)
	#print os.path.basename(snLi)[4:-5], len(lin1), aaa
	names.append(os.path.basename(snLi)[6:-5])
	colN.append(len(lin1))
	aexps.append(aaa)
	
out_name = os.path.join(os.environ["MD10"], 'hlists_MD_1.0Gpc.fits')

redshift = 1./n.array(aexps)-1.
dCom = cosmoMD.comoving_distance(redshift)
ids = n.argsort(redshift)
col0 = fits.Column(name='snap_name'	,format='D', array = n.array(names)[ids] )
col1 = fits.Column(name='N_columns'	,format='I', array = n.array(colN )[ids] )
col2 = fits.Column(name='aexp'		,format='D', array = n.array(aexps)[ids] )
col3 = fits.Column(name='redshift'	,format='D', array = redshift[ids] )
col4 = fits.Column(name='comoving_distance',format='D', array= dCom.value[ids] )
array = 10**9* cosmoMD.age(redshift).value
col5 = fits.Column(name='age_yr'		,format='D', array = array[ids] )
array = cosmoMD.arcsec_per_kpc_comoving(redshift).value*3.6
col6 = fits.Column(name='deg_per_Mpc_comoving', format='D', array = array[ids] )
col7 = fits.Column(name='rho_crit', format='D', array = cosmoMD.critical_density(redshift).to(uu.solMass/(uu.kpc**3)) )
col8 = fits.Column(name='delta_vir', format='D', array = DeltaVir_bn98(redshift) )

#define the table hdu 
hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6, col7, col8])#, col9, col10])
tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
#define the header
prihdr = fits.Header()
prihdr['sim'] = 'MDPL2'
prihdr['author'] = 'JC'
prihdu = fits.PrimaryHDU(header=prihdr)
#writes the file
thdulist = fits.HDUList([prihdu, tb_hdu])
os.system("rm "+out_name)
thdulist.writeto(out_name)

