#  cd pySU/pyMultidark/trunk/bin/fortranfile-0.2.1/

import numpy as n
import os
from os.path import join
from astropy.io import fits
import time
import matplotlib.pyplot as p
box_path=join("..", "MD_1Gpc", "density_field", "Box_HAM_z0.701838_nbar1.350000e-05_QSO.DF.fits.gz")

def plotResidual(box_path, name="residual-position-qso-z07.png"):
	hd = fits.open(box_path)
	L = 1000.
	grid = 2048.
	dx = L / grid
	xr = hd[1].data['x']%dx
	yr = hd[1].data['y']%dx
	zr = hd[1].data['z']%dx
	db=dx/20.
	bins= n.arange(0,dx+db,db)
	p.hist(xr,bins = bins, label='x', histtype='step')
	p.hist(yr,bins = bins, label='y', histtype='step')
	p.hist(zr,bins = bins, label='z', histtype='step')
	p.legend(loc=3)
	p.xlabel('rest(position / 0.48)')
	p.ylabel('count')
	#p.ylim((500, 15000))
	p.xlim((0-db, dx+db))
	#p.yscale('log')
	p.grid()
	p.savefig(join("..", "MD_1Gpc", "density_field", "plots", name))
	p.clf()


box_path=join("..", "MD_1Gpc", "density_field", "Box_HAM_z0.701838_nbar1.350000e-05_QSO.DF.fits.gz")
plotResidual(box_path, name="residual-position-qso-z07.png")

box_path=join("..", "MD_1Gpc", "density_field", "Box_HAM_z0.701838_nbar2.400000e-04_ELG.DF.fits.gz")
plotResidual(box_path, name="residual-position-elg-z07.png")

box_path=join("..", "MD_1Gpc", "density_field", "Box_HAM_z0.701838_nbar1.000000e-04_LRG.DF.fits.gz")
plotResidual(box_path, name="residual-position-lrg-z07.png")

