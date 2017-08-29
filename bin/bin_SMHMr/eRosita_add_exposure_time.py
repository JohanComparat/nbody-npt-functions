import astropy.io.fits as fits
import os
import astropy.wcs as wcs
import numpy as n
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5, Galactic, ICRS, BarycentricTrueEcliptic, HeliocentricTrueEcliptic
import sys

path_2_exposure_map = os.path.join(os.environ['DARKSIM_DIR'], 'observations', "eRosita", "exposure_maps", "eRASS_PeclM27_gal.fits")
texp = fits.open(path_2_exposure_map)
head_wcs = wcs.WCS(texp[0].header)

path_2_h1nh = os.path.join(os.environ['DARKSIM_DIR'], 'observations','h1_nh_LAB.fits')
h1nh = fits.open(path_2_h1nh)
head_wcs_h1nh = wcs.WCS(h1nh[0].header)

path_2_ebv = os.path.join(os.environ['DARKSIM_DIR'], 'observations', 'dust_IREbv.fits')
ebv = fits.open(path_2_ebv)
head_wcs_ebv = wcs.WCS(ebv[0].header)

nhtot = lambda nh1, ebv : nh1 + 2.* 7.2*10**(20) * (1. - n.e**( -nh1 * ebv / (3.*10**(20) ) ) )**(1.1)

path_2_catalog = os.path.join(os.environ['MD10'], "erosita_4most","MDPL2_FluxProj000_ClustersCombinedModel_withHeader.fits")
hd = fits.open(path_2_catalog)
data = hd[1].data

path_2_output = os.path.join(os.environ['MD10'], "erosita_4most","MDPL2_FluxProj000_ClustersCombinedModel_withHeader_withExpTime.fits")

glon = data['galactic_longitude_deg']-texp[0].header['CRPIX1']+0.5 # texp[0].header['CDELT2']/2.
glat = data['galactic_latitude_deg']-texp[0].header['CRPIX2']+0.5 # texp[0].header['CDELT2']/2.
print("ra dec conversion")

c = SkyCoord(glon, glat, frame = Galactic, unit="deg")
out = c.fk5
ecl = c.barycentrictrueecliptic

print("exposure map calculation")
pixX,pixY=head_wcs.wcs_world2pix(glon, glat, 0)
pX_val = pixX.astype('int')
pY_val = pixY.astype('int')
texp_catalog = texp[0].data.T[pX_val,pY_val]

print("ebv calculation")
pixX,pixY=head_wcs_ebv.wcs_world2pix(glon, glat, 0)
pX_val = pixX.astype('int')
pY_val = pixY.astype('int')
ebv_catalog = ebv[0].data.T[pX_val,pY_val]

print("H1nH calculation")
pixX,pixY=head_wcs_h1nh.wcs_world2pix(glon, glat, 0)
pX_val = pixX.astype('int')
pY_val = pixY.astype('int')
h1nh_catalog = h1nh[0].data.T[pX_val,pY_val]

nhtot_catalog = nhtot(h1nh_catalog, ebv_catalog)

exposure_time_col = fits.Column(name='exposure_time', format='D', array=texp_catalog)
number_count_col = fits.Column(name='number_counts', format='D', array=texp_catalog*data['count_rate'])
ebv_col = fits.Column(name='EBV_SFD89', format='D', array=ebv_catalog)
h1nh_col = fits.Column(name='H1NH_Ka05', format='D', array=h1nh_catalog)
nhtot_col = fits.Column(name='nH_total_Wi13', format='D', array=nhtot_catalog)

split_DE_RU_col = fits.Column(name='split', format='L', array = ecl.lat>0 )

ra_col = fits.Column(name='ra', format='D', array=out.ra)
dec_col = fits.Column(name='dec', format='D', array=out.dec)

new_cols  = fits.ColDefs([ exposure_time_col, number_count_col, ra_col, dec_col, ebv_col, h1nh_col, nhtot_col, split_DE_RU_col ])

tb_hdu = fits.BinTableHDU.from_columns( hd[1].columns + new_cols )
prihdr = fits.Header()
prihdr['author'] = 'JC'
prihdu = fits.PrimaryHDU(header=prihdr)
thdulist = fits.HDUList([prihdu, tb_hdu])
os.system("rm "+path_2_output)
thdulist.writeto(path_2_output)


