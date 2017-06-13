import fileinput
import sys
import glob
import os
import numpy as n
import astropy.io.fits as fits

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoNGC = FlatLambdaCDM(H0=68.0*u.km/u.s/u.Mpc, Om0=0.310000, Ob0=0.048)

snList= n.array(glob.glob(os.path.join(os.environ["NUGC_DIR"], "snapshots", "n2gc-*.rockstar")))
snList.sort()

def get_first(path_2_snapshot):
	fl = fileinput.input(path_2_snapshot)
	for line in fl:
		if line[0] == "#" :
			outt = 1 #print line
			if line[1] =="a":
				aexp = float(line[5:])
			if line[1] =="P":
				mass_particle = float(line.split()[-2])
			if line[1] =="B":
				box_length = float(line.split()[-2])
				
		else :
			fl.close()
			return line.split(), aexp, mass_particle, box_length

print "p snapshots"
names, colN, aexps, mps, Lboxs = [], [], [], [], []
for snLi in snList:
	lin1, aaa, mass_particle, box_length = get_first(snLi)
	print aaa, mass_particle, box_length
	print os.path.basename(snLi), len(os.path.basename(snLi))
	names.append(os.path.basename(snLi))
	colN.append(len(lin1))
	aexps.append(aaa)
	mps.append(mass_particle)
	Lboxs.append(box_length)
	
out_name = os.path.join(os.environ["NUGC_DIR"], 'output_NUGC.fits')

redshift = 1./n.array(aexps)-1.
dCom = cosmoNGC.comoving_distance(redshift)
ids = n.argsort(redshift)

col0 = fits.Column(name='snap_name'	,format='A20', array = n.array(names)[ids] )
col1 = fits.Column(name='N_columns'	,format='I', array = n.array(colN )[ids] )
col2 = fits.Column(name='aexp'		,format='D', array = n.array(aexps)[ids] )
col3 = fits.Column(name='redshift'	,format='D', array = redshift[ids] )
col4 = fits.Column(name='comoving_distance',format='D', array= dCom.value[ids] )
array = 10**9* cosmoNGC.age(redshift).value
col5 = fits.Column(name='age_yr'		,format='D', array = array[ids] )
array = cosmoNGC.arcsec_per_kpc_comoving(redshift).value*3.6
col6 = fits.Column(name='deg_per_Mpc_comoving', format='D', array = array[ids] )
col7 = fits.Column(name='box_length', format='D', array = n.array(Lboxs)[ids] )
col8 = fits.Column(name='mass_particle', format='D', array = n.array(mps)[ids] )


#define the table hdu 
hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6, col7, col8])#, col9, col10])
tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
#define the header
prihdr = fits.Header()
prihdr['sim'] = 'NUGC'
prihdr['author'] = 'JC'
prihdu = fits.PrimaryHDU(header=prihdr)
#writes the file
thdulist = fits.HDUList([prihdu, tb_hdu])
os.system("rm "+out_name)
thdulist.writeto(out_name)

