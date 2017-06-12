import numpy as n
zz = 0.
log_Mh = 12.5

# set of parameters / equations
# SMHMr

# SMHMr: exponential cut-off of the SMHMr
nu = lambda z : n.e**(-4./(1+z)**2)

# SMHMr: characteristic halo mass
log10_M1 = lambda z : 11.514 + ( 1.793 * z/(1+z) - 0.251 * z )*nu(z)

# SMHMr: characteristic stellar mass to halo mass ratio
log10_epsilon = lambda z : -1.777 + 0.006 * z/(1+z) * nu(z) + 0.119 * z/(1+z)

# SMHMr: faint end slope
alpha = lambda z : -1.412 - 0.731 * z/(1+z) * nu(z)

# SMHMr: strength of the subpower-law at massive end
delta = lambda z : 3.508 +(- 2.608 * z/(1+z) - 0.043 * z )* nu(z)

# SMHMr: Index of subpower law at massive end
gamma = lambda z : 0.316 +(- 1.319 * z/(1+z) + 0.279 * z )* nu(z)

#SMHMr: Scatter in dex of true stellar mass at fixed halo mass
xi = lambda z : 0.218 +  0.023 * z/(1+z) 

# SMHMr : Systematic offset in stellar masses for active galaxies
kappa = lambda z : 0.045 + 0.155 * z/(1+z)

# SMHMr : Scatter in measured stellar mass at fixed true stellar mass
sigma = lambda z : 0.07 + 0.061 * (z-0.1)

# SFR

# SFR: Characteristic halo mass at which half of stellar mass growth is due to mergers
log10_Mh_ICL = lambda z : 12.515 +  2.503 * z/(1+z) 

# SFR: Systematic offset in stellar masses and SFRs
mu = lambda z : -0.020 - 0.81 * z/(1+z)

# SFR: exponent of the SFR fraction law
beta = lambda z : n.log10(99.)/(16.-log10_Mh_ICL(z))

# SFR: Correlation of SFR to stellar mass at fixed halo mass at a = 0.5
rho_05 = 0.799

# completeness correction as a function of redshift
c_i_z = lambda z : 0.273 / ( 1 + n.e**(1.077 -z))
completeness = lambda x : n.piecewise(x, [x<1, x>=1], [lambda x : 1, lambda x : c_i_z(x) + 1 - c_i_z(1.)]) 

# Fraction of incompleteness due to burstiness (as opposed to dustiness)
b = 0.823

# Model
# Model equations to deduce stellar mass from halo masses

f_x = lambda x, z : delta(z)*(n.log10( 1 + n.e**x))**gamma(z)/(1+n.e**(10**(-x))) - n.log10(10**(alpha(z)*x) + 1)

# Median stellar mass as a function of halo mass and redshift
log10_Mstar = lambda log10_Mhalo, z : log10_epsilon(z) + log10_M1(z) + f_x(log10_Mhalo - log10_M1(z), z) - f_x(0, z)

log10_Mstar_measured = lambda log10_Mstar, z : log10_Mstar + mu(z)
log10_Mstar_measured_active = lambda log10_Mstar, z : log10_Mstar + mu(z) + kappa(z)
f_passive = lambda log10_Mstar, z : 1. / ( 1 + (log10_Mstar_measured(log10_Mstar, z)/10**(10.2+0.5*z))**(-1.3) )
f_loss = lambda time : 0.05*n.log(1+time/1.4) # time in Myr
# star formation rate fraction
f_SFR = lambda log10_Mhalo, z: 1/(1 + (10**(log10_Mhalo-log10_Mh_ICL(z)))**beta(z))
