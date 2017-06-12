# overall python packages
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

# dedicated packages
import ClusterScalingRelations
cl = ClusterScalingRelations.ClusterScalingRelations_Mantz2016()
import StellarMass
import XrayLuminosity
xr = XrayLuminosity.XrayLuminosity()

import matplotlib.pyplot as p


smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
ll_dir = os.path.join(os.environ['DATA_DIR'], 'spm', 'literature')
path_ilbert13_SMF = os.path.join(ll_dir, "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(os.path.join( ll_dir, "ilbert_2013_mass_function_params.txt"), unpack=True)

smf01 = lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
print 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0]
smf08 = lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )
print 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] 


logMs = n.arange(6.5,12.5,0.01)
mbins = n.arange(8,12.5,0.25)

p.figure(1, (6,6))
    
p.plot(mbins, n.log10(smf01(10**mbins)), label='Ilbert 13, 0.2<z<0.5', ls='dashed')
p.plot(mbins, n.log10(smf08(10**mbins)), label='Ilbert 13, 0.8<z<1.1', ls='dashed')
    
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 0.3)) for logMs_i in logMs]) , label='z=0.3')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 0.9)) for logMs_i in logMs]) , label='z=0.9')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 1.2)) for logMs_i in logMs]) , label='z=1.2')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 1.8)) for logMs_i in logMs]) , label='z=1.8')
p.plot(logMs, n.array([n.log10(xr.Phi_stellar_mass(logMs_i, 2.5)) for logMs_i in logMs]) , label='z=2.5')
p.xlabel('stellar mass')
p.ylabel('log Phi stellar mass')
p.title('AGN HGMF')
p.xlim((9.5, 12.5))
p.ylim((-6,-2))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig('/home/comparat/data/eRoMok/BO12_fig3.png')
p.clf()


sys.exit()
print " set up box, and redshift "
#MD 1 hlist_0.74980_SAM_Nb_0.fits
#MD 25 hlist_0.75440_SAM_Nb_10.fits

z = 1./0.74230 -1.
# set up the stellar mass computation
sm = StellarMass.StellarMass()
mhs = n.logspace(10,15,100)
Nhalo = len(mhs)

ratio = sm.SMHMr(mhs,0.)
stellar_mass = sm.meanSM(mhs,0.)

# set up the x ray lambda SAR
logMs = n.arange(6.5,12.5,0.01)
cdfs_interpolations = []
cdfs_interpolations_maxs = []
XXS = n.arange(32,36.1,0.1)
for jj,mass in enumerate(logMs):
    pd = lambda ll : xr.psi_log(ll, logM=mass, z=z)
    norming = quad( pd, 32, 36)[0]
    cdfs_interpolations.append( interp1d(n.array([quad( pd, 32, X)[0] for X in XXS ])/norming, XXS) )
    cdfs_interpolations_maxs.append( norming*50. )

cdfs_interpolations = n.array(cdfs_interpolations)
cdfs_interpolations_maxs = n.array(cdfs_interpolations_maxs)

Mgal_mvir_Mo13 = norm.rvs( loc = sm.meanSM(mhs, z), scale = 0.15 )
randomX = n.random.rand(len(Mgal_mvir_Mo13))
indexes = n.searchsorted(logMs,Mgal_mvir_Mo13)
lambda_sar_Bo16 = n.array([ cdfs_interpolations[indexes[ii]](randomX[ii]) for ii in range(Nhalo) ])
active_gn = n.array([ cdfs_interpolations_maxs[indexes[ii]] > randomX[ii] for ii in range(Nhalo) ])

Mgal_m200c_Mo13 = norm.rvs( loc = sm.meanSM(mhs, z), scale = 0.15 )
randomY = n.random.rand(len(Mgal_m200c_Mo13))
indexesY = n.searchsorted(logMs,Mgal_m200c_Mo13)
lambda_sar_Bo16_m200c = n.array([ cdfs_interpolations[indexesY[ii]](randomY[ii]) for ii in range(Nhalo) ])
active_gn_m200c = n.array([ cdfs_interpolations_maxs[indexesY[ii]] > randomY[ii] for ii in range(Nhalo) ])

# columns related to clusters
Mgas_cluster   = n.log10(cl.logM500_to_logMgas(n.log10(mhs), z))
kT_cluster     = cl.logM500_to_kT(n.log10(mhs), z)
Lx_bol_cluster = n.log10(cl.logM500_to_L(n.log10(mhs), z))
Lx_ce_cluster  = n.log10(cl.logM500_to_Lce(n.log10(mhs), z))

logLSAR =n.arange(32,35,0.1)
z=0.5
p.figure(1, (6,6))
p.plot(logLSAR, n.log10(xr.f_lambda_sar(9.5, z, logLSAR) ), label='9.5')
p.plot(logLSAR, n.log10(xr.f_lambda_sar(10., z, logLSAR) ), label='10.')
p.plot(logLSAR, n.log10(xr.f_lambda_sar(10.5, z, logLSAR)), label='10.5')
p.plot(logLSAR, n.log10(xr.f_lambda_sar(11., z, logLSAR) ), label='11.')
p.plot(logLSAR, n.log10(xr.f_lambda_sar(11.5, z, logLSAR)), label='11.5')
p.plot(logLSAR, n.log10(xr.f_lambda_sar(12., z, logLSAR) ), label='12')
p.xlabel('lambda SAR')
p.ylabel('f lambda SAR')
p.title('z=0.3')
p.xlim((32,35.5))
p.ylim((-5,4))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig('/home/comparat/data/eRoMok/BO12_fig4_up.png')
p.clf()

p.figure(1, (6,6))
p.plot(logLSAR, n.array([n.log10(xr.Phi_lambda_SAR(logLSAR_i, 0.3)) for logLSAR_i in logLSAR]) , label='z=0.3')
p.plot(logLSAR, n.array([n.log10(xr.Phi_lambda_SAR(logLSAR_i, 0.6)) for logLSAR_i in logLSAR]) , label='z=0.6')
p.plot(logLSAR, n.array([n.log10(xr.Phi_lambda_SAR(logLSAR_i, 0.9)) for logLSAR_i in logLSAR]) , label='z=0.9')
p.plot(logLSAR, n.array([n.log10(xr.Phi_lambda_SAR(logLSAR_i, 1.2)) for logLSAR_i in logLSAR]) , label='z=1.2')
p.plot(logLSAR, n.array([n.log10(xr.Phi_lambda_SAR(logLSAR_i, 1.5)) for logLSAR_i in logLSAR]) , label='z=1.5')
p.plot(logLSAR, n.array([n.log10(xr.Phi_lambda_SAR(logLSAR_i, 1.8)) for logLSAR_i in logLSAR]) , label='z=1.8')
p.plot(logLSAR, n.array([n.log10(xr.Phi_lambda_SAR(logLSAR_i, 2.2)) for logLSAR_i in logLSAR]) , label='z=2.2')
p.plot(logLSAR, n.array([n.log10(xr.Phi_lambda_SAR(logLSAR_i, 2.5)) for logLSAR_i in logLSAR]) , label='z=2.5')
p.xlabel('lambda SAR')
p.ylabel('Phi lambda SAR')
p.title('SARDF')
p.xlim((32,35.5))
p.ylim((-8,-3))
p.grid()
p.legend(loc=0, frameon=False)
p.savefig('/home/comparat/data/eRoMok/BO12_fig4_low.png')
p.clf()
