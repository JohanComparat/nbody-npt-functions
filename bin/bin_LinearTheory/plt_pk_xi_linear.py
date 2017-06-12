import astropy.units as uu
import astropy.cosmology as co
aa = co.Planck13
from scipy.interpolate import interp1d
import numpy as n
import matplotlib
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p
import glob
import sys
from os.path import join
from scipy.optimize import minimize

zmin = 0.
zmax = 3.5
NDecimal = 3

# defining directories :
dir = ".." #join("D:","data","MultiDark")
zList_file =  join(dir, "z-list-all-boxes.txt") 
z0 = n.loadtxt(zList_file,unpack=True)
pklist = n.array(glob.glob(join(dir, "Pk_DM_CLASS", "MD_z*_xi.dat")))


p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
for j, pkfile in enumerate(pklist):
	print pkfile
	snapNum = pkfile.split('_')[-2][1:]
	zSnap = z0[ int(snapNum) - 1 ]
	rr, xi = n.loadtxt(pkfile, unpack=True)
	if zSnap>=zmin and zSnap<=zmax :
		sc1=p.scatter(rr, rr*rr*xi, c=n.log10(1+zSnap*n.ones_like(rr)),s=5, marker='o', rasterized=True, vmin = n.log10(1+zmin), vmax =n.log10(1+ zmax))
		sc1.set_edgecolor('face')

cb = p.colorbar(shrink=0.8)
cb.set_label("log(1+z)")
zs = n.hstack((0,0.15,0.3,n.arange(0.5,5,0.5)))
cb.set_ticks(n.log10(1+zs))
cb.set_ticklabels(zs)
p.xlabel(r'$r$')
p.ylabel(r'$r^2\xi(r)$')
p.ylim((-15,50))
p.xlim((0,200))
#p.xscale('log')
#p.yscale('log')
p.grid()

p.savefig(join("..","PK_DM_CLASS", "plots", "XI.png"))
p.clf()

pklist = n.array(glob.glob(join(dir, "Pk_DM_CLASS", "MD_z*_pk.dat")))

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
for j, pkfile in enumerate(pklist):
	print pkfile
	snapNum = pkfile.split('_')[-2][1:]
	zSnap = z0[ int(snapNum) - 1 ]
	ks_lin_0, pks_lin_0 = n.loadtxt(pkfile, unpack=True)
	if zSnap>=zmin and zSnap<=zmax :
		sc1=p.scatter(n.log10(ks_lin_0), n.log10(pks_lin_0), c=n.log10(1+zSnap*n.ones_like(ks_lin_0)),s=5, rasterized=True, vmin = n.log10(1+zmin), vmax =n.log10(1+ zmax))
		sc1.set_edgecolor('face')

cb = p.colorbar(shrink=0.8)
cb.set_label("log(1+z)")
zs = n.hstack((0,0.15,0.3,n.arange(0.5,5,0.5)))
cb.set_ticks(n.log10(1+zs))
cb.set_ticklabels(zs)
p.xlabel(r'$\log(k)$')
p.ylabel(r'$\log(P(k))$')
p.ylim((-2,5))
p.xlim((-3, 1))
#p.xscale('log')
#p.yscale('log')
p.grid()

p.savefig(join("..","PK_DM_CLASS", "plots", "PK.png"))
p.clf()
