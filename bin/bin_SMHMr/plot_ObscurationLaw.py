import numpy as n
from scipy.stats import norm
from scipy.integrate import quad
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
import glob
import astropy.io.fits as fits
import os
import time
import numpy as n
import sys

data=n.loadtxt("/data17s/darksim/software/nbody-npt-functions/data/AGN/fraction_observed_by_erosita_due_2_obscuration.txt")

out_dir = os.path.join(os.path.join(os.environ['MD10'],"results","obscuration_agn", "images"))
	
# compare the stellar mass function measured to the Ilbert function
# take the AGN HGMF model

p.figure(1, (6,6))
p.scatter(n.log10(data.T[0]), data.T[1], c=n.log10(data.T[2]), s=10,edgecolor='face')
cb=p.colorbar()
cb.set_label('log10(unobscured fraction)')
p.xlabel(r'$\log_{10}(1+z)$')
p.ylabel(r'$\log_{10}(Nh [cm^{-2}])$')
p.xlim((-3.01, 0.7))
p.ylim((20,26))
p.title('M. Brightman torus model')
p.grid()
p.savefig(os.path.join(out_dir, "AGN_obscuration.png"))
p.clf()

os.system("cp $MD10/results/obscuration_agn/images/AGN_obscuration.png ~/wwwDir/eRoMok/")

