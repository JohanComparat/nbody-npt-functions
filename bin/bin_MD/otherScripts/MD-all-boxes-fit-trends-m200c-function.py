import astropy.units as uu
import astropy.cosmology as co
aa = co.Planck13
#aah = co.FlatLambdaCDM(H0=100.0 *uu.km / (uu.Mpc *uu.s), Om0=0.307, Tcmb0=2.725 *uu.K, Neff=3.05, m_nu=[ 0.  ,  0. ,   0.06]*uu.eV, Ob0=0.0483)
#rhom = aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aa.critical_density0.to(uu.solMass*uu.Mpc**-3).value
#aah.critical_density0.to(uu.solMass*uu.Mpc**-3).value

from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
import numpy as n
import matplotlib
#matplotlib.use('pdf')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle
from os.path import join
from scipy.optimize import minimize
import glob

#saundersFct=lambda v, A, logv0,alpha,sig : 10**A * ( v /10**logv0)**(alpha) * n.e**( -n.log10( 1 + v /10**logv0)/(2*sig**2.))
#schechterFct=lambda v, A, logv0,alpha, sig : 10**A * ( v /10**logv0)**(alpha) * n.e**( - v / 10**logv0 /(2*sig**2.) )
#ps * (10**logl/10**logls)**(a+1) * n.e**(-10**logl/10**logls)
#doublePL=lambda v,A,logv0 ,a,b: 10**A * 10**( (1+a)*( v - 10**logv0) + b )

mf = lambda v, A, v0, alpha, beta : 10**A * (v/10**v0)**beta * n.e**(- (v/10**v0)**alpha )

# limits at z0
limits_04 = [1e10, 5e12]
limits_10 = [5e11, 5e13]
limits_25 = [5e12, 5e14]
limits_40 = [1e13, 5e15]
zmin = -0.1
zmax = 4

NDecimal = 3

dir = join("D:","\data","BigMD","M200cFunction")

qty_limits = n.array([limits_04, limits_10, limits_25, limits_40])
volume_boxes =  n.array([400.**3., 1000**3., 2500**3., 4000.**3.])

cos = "Central" # centrak or satellite ?
qty = "M200c"

print "we consider the ",type,"of",qty,"of", cos
print "in the redshift range",zmin,zmax
#print zList_files
print "within the following limits for each box",qty_limits
print "each box has a volume of",volume_boxes, "Mpc3/h3"



xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir,"hist-Central-M200c_ALL_cumulative_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir,"hist-Central-M200c_ALL_cumulative_MD_1Gpc.dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir,"hist-Central-M200c_ALL_cumulative_MD_2.5Gpc.dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir,"hist-Central-M200c_ALL_cumulative_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= zmin) & (z_04 <= zmax)
s_10 = (z_10 >= zmin) & (z_10 <= zmax)
s_25 = (z_25 >= zmin) & (z_25 <= zmax)
s_40 = (z_40 >= zmin) & (z_40 <= zmax)

redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)

fitList = n.array(glob.glob(join(dir,"MF-4params-z-*-fit")))

data = n.empty((7,10))

data[0] = n.loadtxt("D:\BigMD\M200cFunction\MF-4params-z-0.0-fit.txt")
data[1] = n.loadtxt("D:\BigMD\M200cFunction\MF-4params-z-0.105225074332-fit.txt")
data[2] = n.loadtxt("D:\BigMD\M200cFunction\MF-4params-z-0.336691790552-fit.txt")
data[3] = n.loadtxt("D:\BigMD\M200cFunction\MF-4params-z-0.673566901501-fit.txt")
data[4] = n.loadtxt("D:\BigMD\M200cFunction\MF-4params-z-1.009782226-fit.txt")
data[5] = n.loadtxt("D:\BigMD\M200cFunction\MF-4params-z-1.46051122015-fit.txt")
#data[5] = n.loadtxt("D:\BigMD\M200cFunction\MF-4params-z-2.10636649588-fit.txt")
data[6] = n.loadtxt("D:\BigMD\M200cFunction\MF-4params-z-2.91314105144-fit.txt")

A0,A0_err,vcut0, vcut0_err, a0, a0_err, b0, b0_err, z, z_err = n.transpose(data)

pl1= lambda x, p0, p1 : p0 + p1*x 
pl2 = lambda x, p0, p1,p2 : p0 + p1*x + p2*x*x


p0 = n.array([-3.9,1])

def fitPL1modelToParameters(x, y, x_err, y_err, p0, ylab):
	model = pl1
	popt, cov = curve_fit(model, x, y, sigma=y_err, p0 = p0 , maxfev = 5000000)
	print popt, cov
	xModel = n.arange(x.min()*0.9,x.max()*1.1,0.1)
	yModel = model(xModel,popt[0],popt[1])
	print xModel, yModel
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.errorbar(x, y, xerr= x_err,yerr = y_err, fmt='none')
	p.plot(xModel, yModel,'k--',label="model")
	p.xlabel('z')
	p.ylabel(ylab) # log$_{10}[ n(>M)]')
	p.title(str(n.round(popt[0],NDecimal))+" "+str(n.round(popt[1],NDecimal)) )
	p.grid()
	p.show()

	
def fitPL2modelToParameters(x, y, x_err, y_err, p0, ylab):
	model = pl2
	popt, cov = curve_fit(model, x, y, sigma=y_err, p0 = p0 , maxfev = 5000000)
	print popt, cov
	xModel = n.arange(x.min()*0.9,x.max()*1.1,0.1)
	yModel = model(xModel,popt[0],popt[1],popt[2])
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.errorbar(x, y, xerr= x_err,yerr = y_err, fmt='none')
	p.plot(xModel, yModel,'k--',label="model")
	p.xlabel('z')
	p.ylabel(ylab) # log$_{10}[ n(>M)]')
	p.title(str(n.round(popt[0],NDecimal))+" "+str(n.round(popt[1],NDecimal))+" "+str(n.round(popt[2],NDecimal)) )
	p.grid()
	p.show()	
	

def fitPL2AncheredmodelToParameters(x, y, x_err, y_err, pGuess,p0, ylab):
	model = lambda x, p1, p2 : pl2(x, p0, p1, p2)
	popt, cov = curve_fit(model, x, y, sigma=y_err, p0 = pGuess , maxfev = 5000000)
	print popt, cov
	xModel = n.arange(x.min()*0.9,x.max()*1.1,0.1)
	yModel = model(xModel,popt[0],popt[1])
	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.errorbar(x, y, xerr= x_err,yerr = y_err, fmt='none')
	p.plot(xModel, yModel,'k--',label="model")
	p.xlabel('z')
	p.ylabel(ylab) # log$_{10}[ n(>M)]')
	p.title(str(n.round(p0,NDecimal))+" "+str(n.round(popt[0],NDecimal))+" "+str(n.round(popt[1],NDecimal)) )
	p.grid()
	p.show()	
		
p0 = n.array([0,0])

fitPL2AncheredmodelToParameters(z, A0, z_err, A0_err, p0, -4.009,r'$A_0$')
fitPL2AncheredmodelToParameters(z, vcut0, z_err, vcut0_err, p0,13.826,r'$M_{cut}$')
fitPL2AncheredmodelToParameters(z, a0, z_err, a0_err, p0,0.578,r'$\alpha$')
fitPL2AncheredmodelToParameters(z, b0, z_err, b0_err, p0,-0.876,r'$\beta$')

sys.exit()

fitPL1modelToParameters(z, A0, z_err, A0_err, p0,r'$A_0$')
fitPL1modelToParameters(z, vcut0, z_err, vcut0_err, p0,r'$M_{cut}$')
fitPL1modelToParameters(z, a0, z_err, a0_err, p0,r'$\alpha$')
fitPL1modelToParameters(z, b0, z_err, b0_err, p0,r'$\beta$')

p0 = n.array([-4,0,0])

fitPL2modelToParameters(z, A0, z_err, A0_err, p0,r'$A_0$')
fitPL2modelToParameters(z, vcut0, z_err, vcut0_err, p0,r'$M_{cut}$')
fitPL2modelToParameters(z, a0, z_err, a0_err, p0,r'$\alpha$')
fitPL2modelToParameters(z, b0, z_err, b0_err, p0,r'$\beta$')

sys.exit()
################################ Plot cumulative halo mass function and model at z=0  ################################

def fitSingleRedshiftM200Function(z0,z1):

	s_04 = (z_04 < z1) & (z_04 > z0 )
	s_10 = (z_10 < z1) & (z_10 > z0 )
	s_25 = (z_25 < z1) & (z_25 > z0 )
	s_40 = (z_40 < z1) & (z_40 > z0 )

	redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
	print "all redshifts available:", set(redshift)
	M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
	print "min and max masses available:", n.min(M200c), n.max(M200c)
	yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
	print "min and max Y available:", n.min(yData), n.max(yData)
	yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
	print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)

	p0 = n.array([-3.5, 12.5, 0.8, -1.])

	vf = lambda v, A, v0, alpha, beta : n.log10( 10**A * (10**v/10**v0)**beta * n.e**(- (10**v/10**v0)**alpha ) )

	# with curve fit
	print "with curve fit"
	popt, cov = curve_fit(vf, M200c, yData,sigma=0.025*n.ones_like(yData), p0 = p0 , maxfev = 5000000)
	print popt, cov
	A0, vcut0, a0, b0 = n.round(popt,NDecimal)

	print "redshift 0 model for the M200c cumulative function :"
	print "A(z=0) & = "+str(A0)+ r"\pm ", n.round(cov[0][0]**0.5, NDecimal), r'\\'
	print r" M_{200c}^{cut}(z=0) & = "+str(vcut0)+ r"\pm ", n.round(cov[1][1]**0.5, NDecimal), r'\\'
	print r" \alpha(z=0) & = "+str(a0)+ r"\pm ", n.round(cov[2][2]**0.5, NDecimal), r'\\'
	print r" \beta(z=0) & = "+str(b0)+ r"\pm ", n.round(cov[3][3]**0.5, NDecimal), r'\\'

	f=open("MF-4params-z-"+str(n.mean(redshift))+"-fit.txt",'w')
	n.savetxt(f,n.array([A0,cov[0][0]**0.5,vcut0,cov[1][1]**0.5, a0, cov[2][2]**0.5, b0, cov[3][3]**0.5, n.mean(redshift), n.std(redshift)]))
	f.close()

	p.figure(0,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.plot(n.log10(xData_04[s_04][::3]), n.log10(yData_04[s_04][::3]), marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)
	p.plot(n.log10(xData_10[s_10][::3]), n.log10(yData_10[s_10][::3]), marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)
	p.plot(n.log10(xData_25[s_25][::3]), n.log10(yData_25[s_25][::3]), marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)
	p.plot(n.log10(xData_40[s_40][::3]), n.log10(yData_40[s_40][::3]), marker ='p', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)
	xModel = n.arange(n.min(M200c),n.max(M200c),0.1)
	yModel_CF = vf(xModel, A0, vcut0, a0, b0)
	p.plot(xModel, yModel_CF,'k--',label="model")
	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r' n(>M)') # log$_{10}[ n(>M)]')
	p.legend(loc=3)
	p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0))
	p.grid()

	yModel_CF = vf(M200c, A0, vcut0, a0, b0)

	p.figure(1,(6,6))
	p.axes([0.17,0.17,0.75,0.75])
	p.plot(M200c,10**yData/10**yModel_CF,'bo')
	p.axhline(1)

	p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
	p.ylabel(r' n(>M) data / model') # log$_{10}[ n(>M)]')
	p.legend(loc=3)
	p.ylim((.95,1.05))
	p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0))
	p.grid()
	p.show()


fitSingleRedshiftM200Function(-0.1,0.01)
fitSingleRedshiftM200Function(0.091,0.15)
fitSingleRedshiftM200Function(0.25,0.35)
fitSingleRedshiftM200Function(0.65,0.7)
fitSingleRedshiftM200Function(0.99,1.1)
fitSingleRedshiftM200Function(1.4,1.5)
fitSingleRedshiftM200Function(2.,2.2)
fitSingleRedshiftM200Function(2.8,3.)



sys.exit()






















"""
p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])

sc1=p.scatter(M200c, yData, c=redshift,s=10, marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r' n(>M)') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.grid()
p.show()
p.savefig("M200c-cumulative-function-data.pdf")
p.clf()
"""
################ now does the fits at redshift 0

s_04 = (z_04 < 0.01) #& (z_04 > 2. )
s_10 = (z_10 < 0.01) #& (z_10 > 2. )
s_25 = (z_25 < 0.01) #& (z_25 > 2. )
s_40 = (z_40 < 0.01) #& (z_40 > 2. )

redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)


p0 = n.array([-3.5, 12.5, 0.8, -1.])

vf = lambda v, A, v0, alpha, beta : n.log10( 10**A * (10**v/10**v0)**beta * n.e**(- (10**v/10**v0)**alpha ) )

# with curve fit
print "with curve fit"
popt, cov = curve_fit(vf, M200c, yData,sigma=0.025*n.ones_like(yData), p0 = p0 , maxfev = 5000000)
print popt, cov
A0, vcut0, a0, b0 = n.round(popt,NDecimal)

print "redshift 0 model for the M200c cumulative function :"
print "A(z=0) & = "+str(A0)+ r"\pm ", n.round(cov[0][0]**0.5, NDecimal), r'\\'
print r" M_{200c}^{cut}(z=0) & = "+str(vcut0)+ r"\pm ", n.round(cov[1][1]**0.5, NDecimal), r'\\'
print r" \alpha(z=0) & = "+str(a0)+ r"\pm ", n.round(cov[2][2]**0.5, NDecimal), r'\\'
print r" \beta(z=0) & = "+str(b0)+ r"\pm ", n.round(cov[3][3]**0.5, NDecimal), r'\\'


f=open("MF-4params-z-"+str(n.mean(redshift))+"-fit.txt",'w')
n.savetxt(f,n.array([A0,cov[0][0]**0.5,vcut0,cov[1][1]**0.5, a0, cov[2][2]**0.5, b0, cov[3][3]**0.5, n.mean(redshift), n.std(redshift)]))
f.close()

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.plot(n.log10(xData_04[s_04][::3]), n.log10(yData_04[s_04][::3]), marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)
p.plot(n.log10(xData_10[s_10][::3]), n.log10(yData_10[s_10][::3]), marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)
p.plot(n.log10(xData_25[s_25][::3]), n.log10(yData_25[s_25][::3]), marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)
p.plot(n.log10(xData_40[s_40][::3]), n.log10(yData_40[s_40][::3]), marker ='p', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)
xModel = n.arange(n.min(M200c),n.max(M200c),0.1)
yModel_CF = vf(xModel, A0, vcut0, a0, b0)
p.plot(xModel, yModel_CF,'k--',label="model")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r' n(>M)') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0))
p.grid()

yModel_CF = vf(M200c, A0, vcut0, a0, b0)

p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.plot(M200c,10**yData/10**yModel_CF,'bo')
p.axhline(1)

p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r' n(>M) data / model') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.ylim((.95,1.05))
p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0))
p.grid()
p.show()


################ now does the fits at redshift 2.

s_04 = (z_04 < 2.8) & (z_04 > 2. )
s_10 = (z_10 < 2.8) & (z_10 > 2. )
s_25 = (z_25 < 2.8) & (z_25 > 2. )
s_40 = (z_40 < 2.8) & (z_40 > 2. )

redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)
redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)


p0 = n.array([-3.5, 12.5, 0.8, -1.])

vf = lambda v, A, v0, alpha, beta : n.log10( 10**A * (10**v/10**v0)**beta * n.e**(- (10**v/10**v0)**alpha ) )

# with curve fit
print "with curve fit"
popt, cov = curve_fit(vf, M200c, yData,sigma=0.025*n.ones_like(yData), p0 = p0 , maxfev = 5000000)
print popt, cov
A0, vcut0, a0, b0 = n.round(popt,NDecimal)

print "redshift 0 model for the M200c cumulative function :"
print "A(z=0) & = "+str(A0)+ r"\pm ", n.round(cov[0][0]**0.5, NDecimal), r'\\'
print r" M_{200c}^{cut}(z=0) & = "+str(vcut0)+ r"\pm ", n.round(cov[1][1]**0.5, NDecimal), r'\\'
print r" \alpha(z=0) & = "+str(a0)+ r"\pm ", n.round(cov[2][2]**0.5, NDecimal), r'\\'
print r" \beta(z=0) & = "+str(b0)+ r"\pm ", n.round(cov[3][3]**0.5, NDecimal), r'\\'


f=open("MF-4params-z-"+str(n.mean(redshift))+"-fit.txt",'w')
n.savetxt(f,n.array([A0,cov[0][0]**0.5,vcut0,cov[1][1]**0.5, a0, cov[2][2]**0.5, b0, cov[3][3]**0.5, n.mean(redshift), n.std(redshift)]))
f.close()

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.plot(n.log10(xData_04[s_04][::3]), n.log10(yData_04[s_04][::3]), marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)
p.plot(n.log10(xData_10[s_10][::3]), n.log10(yData_10[s_10][::3]), marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)
p.plot(n.log10(xData_25[s_25][::3]), n.log10(yData_25[s_25][::3]), marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)
p.plot(n.log10(xData_40[s_40][::3]), n.log10(yData_40[s_40][::3]), marker ='p', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)
xModel = n.arange(n.min(M200c),n.max(M200c),0.1)
yModel_CF = vf(xModel, A0, vcut0, a0, b0)
p.plot(xModel, yModel_CF,'k--',label="model")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r' n(>M)') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0))
p.grid()
p.show()

################ now does the fits at redshift 2.8

s_04 = (z_04 >= 2.8) #& (z_04 > 2. )
s_10 = (z_10 >= 2.8) #& (z_10 > 2. )
s_25 = (z_25 >= 2.8) #& (z_25 > 2. )
s_40 = (z_40 >= 2.8) #& (z_40 > 2. )

redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)
redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)


p0 = n.array([-3.5, 12.5, 0.8, -1.])

vf = lambda v, A, v0, alpha, beta : n.log10( 10**A * (10**v/10**v0)**beta * n.e**(- (10**v/10**v0)**alpha ) )

# with curve fit
print "with curve fit"
popt, cov = curve_fit(vf, M200c, yData,sigma=0.025*n.ones_like(yData), p0 = p0 , maxfev = 5000000)
print popt, cov
A0, vcut0, a0, b0 = n.round(popt,NDecimal)

print "redshift 0 model for the M200c cumulative function :"
print "A(z=0) & = "+str(A0)+ r"\pm ", n.round(cov[0][0]**0.5, NDecimal), r'\\'
print r" M_{200c}^{cut}(z=0) & = "+str(vcut0)+ r"\pm ", n.round(cov[1][1]**0.5, NDecimal), r'\\'
print r" \alpha(z=0) & = "+str(a0)+ r"\pm ", n.round(cov[2][2]**0.5, NDecimal), r'\\'
print r" \beta(z=0) & = "+str(b0)+ r"\pm ", n.round(cov[3][3]**0.5, NDecimal), r'\\'


f=open("MF-4params-z-"+str(n.mean(redshift))+"-fit.txt",'w')
n.savetxt(f,n.array([A0,cov[0][0]**0.5,vcut0,cov[1][1]**0.5, a0, cov[2][2]**0.5, b0, cov[3][3]**0.5, n.mean(redshift), n.std(redshift)]))
f.close()

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])
p.plot(n.log10(xData_04[s_04][::3]), n.log10(yData_04[s_04][::3]), marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)
p.plot(n.log10(xData_10[s_10][::3]), n.log10(yData_10[s_10][::3]), marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)
p.plot(n.log10(xData_25[s_25][::3]), n.log10(yData_25[s_25][::3]), marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)
p.plot(n.log10(xData_40[s_40][::3]), n.log10(yData_40[s_40][::3]), marker ='p', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)
xModel = n.arange(n.min(M200c),n.max(M200c),0.1)
yModel_CF = vf(xModel, A0, vcut0, a0, b0)
p.plot(xModel, yModel_CF,'k--',label="model")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r' n(>M)') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.title(str(n.round(n.mean(redshift),NDecimal))+" "+str(A0)+" "+str(vcut0)+" "+str(a0)+" "+str(b0))
p.grid()
p.show()



sys.exit()

























################################ Model Fits on the cumulative function, evolution with redshift ################################

xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir,"hist-Central-M200c_ALL_cumulative_MD_0.4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir,"hist-Central-M200c_ALL_cumulative_MD_1Gpc.dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir,"hist-Central-M200c_ALL_cumulative_MD_2.5Gpc.dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir,"hist-Central-M200c_ALL_cumulative_MD_4Gpc.dat"),unpack=True)

s_04 = (z_04 >= 0.0) & (z_04 <= zmax)
s_10 = (z_10 >= 0.0) & (z_10 <= zmax)
s_25 = (z_25 >= 0.0) & (z_25 <= zmax)
s_40 = (z_40 >= 0.0) & (z_40 <= zmax)

redshift = n.hstack(( z_04[s_04], z_10[s_10], z_25[s_25], z_40[s_40]))
print "all redshifts available:", set(redshift)
M200c = n.log10(n.hstack(( xData_04[s_04], xData_10[s_10], xData_25[s_25], xData_40[s_40])))
print "min and max masses available:", n.min(M200c), n.max(M200c)
yData = n.log10(n.hstack(( yData_04[s_04], yData_10[s_10], yData_25[s_25], yData_40[s_40])))
print "min and max Y available:", n.min(yData), n.max(yData)
yDataErr = abs(n.hstack(( yDataErr_04[s_04], yDataErr_10[s_10], yDataErr_25[s_25], yDataErr_40[s_40])) / yData)
print "min and max Y error available:", n.min(yDataErr), n.max(yDataErr)

#vfG = lambda v, z, A0, A1, vcut0, vcut1, b0, b1 : 10**(A0 + A1 * z) * (1+ (v/10**(vcut0 + vcut1 * z))**(b0 + b1 * z) )* n.e**(- (v/10**(vcut0 + vcut1 * z))**(a0 ) )
vfG = lambda v, z, A1, vcut1, a1, b1 : n.log10( 10**(A0 + A1 * z) * (1+ (10**v/10**(vcut0 + vcut1 * z))**(b0 + b1 * z) )* n.e**(- (10**v/10**(vcut0 + vcut1 * z))**(a0 +a1*z) ))
vfGbis = lambda v, z, ps : vfG(v,z,ps[0],ps[1],ps[2],ps[3])
chi2fun = lambda ps : n.sum( (vfGbis(M200c,redshift,ps) - yData)**2. / (n.ones_like(yData)*0.05)**2. )/len(yData)
p1 = n.array([ 0., 0., 0., 0.])
print "looks for the optimum parameters"
res = minimize(chi2fun, p1, method='Powell',options={'xtol': 1e-6, 'disp': True, 'maxiter' : 50000000, 'nfev': 1800000})

print "ndof=",len(yData)
print res
A1, vcut1, a1, b1 = n.round(res.x,4)

print "A(z) & = "+str(A0)+" + "+str(A1)+r'\times z \\'
print r" M_{cut}(z) & = "+str(vcut0)+" + "+str(vcut1)+r'\times z \\'
print r" \alpha(z) & = "+str(a0)+" + "+str(a1)+r'\times z \\' #+ '+str(a2)+r'\times z^2 \\'
print r" \beta(z) & = "+str(b0)+" + "+str(b1)+r'\times z \\'


vfG = lambda v, z, vcut1, b1 : n.log10( 10**(A0 ) * (1+ (10**v/10**(vcut0 + vcut1 * z))**(b0 + b1 * z) )* n.e**(- (10**v/10**(vcut0 + vcut1 * z))**(a0 ) ))
vfGbis = lambda v, z, ps : vfG(v,z,ps[0],ps[1])
chi2fun = lambda ps : n.sum( (vfGbis(M200c,redshift,ps) - yData)**2. / (n.ones_like(yData)*0.05)**2. )/len(yData)

p1 = n.array([ 0.,  0.])

print "looks for the optimum parameters"
res = minimize(chi2fun, p1, method='Powell',options={'xtol': 1e-6, 'disp': True, 'maxiter' : 50000000, 'nfev': 1800000})

print "ndof=",len(yData)
print res
A1, vcut1, a1, b1 = n.round(res.x,4)

print "A(z) & = "+str(A0)+" + "+str(A1)+r'\times z \\'
print r" M_{cut}(z) & = "+str(vcut0)+" + "+str(vcut1)+r'\times z \\'
print r" \alpha(z) & = "+str(a0)+" + "+str(a1)+r'\times z \\' #+ '+str(a2)+r'\times z^2 \\'
print r" \beta(z) & = "+str(b0)+" + "+str(b1)+r'\times z \\'

# now outputs the model
xModel = n.arange(n.min(M200c),15,0.1)

X,Y = n.meshgrid(xModel,n.arange(0,n.max(redshift)+0.025,0.025))

Z = vfGbis(X,Y,res.x)

n.savetxt(join(dir,"M200c-cumulative-function-best_fit.txt"),n.transpose([n.hstack((X)), n.hstack((Y)), n.hstack((Z))]) )

#######################################################
# now plots the results of the fit
print "now plots the results of the fit"

vmax_mod, z_mod, n_mod = n.loadtxt(join(dir,"M200c-cumulative-function-best_fit.txt"), unpack=True)


p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])

sc1=p.scatter(vmax_mod, n_mod, c=z_mod,s=10, marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
cb = p.colorbar(shrink=0.8)
cb.set_label("redshift")
p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r' n(>M)') # log$_{10}[ n(>M)]')
p.legend(loc=3)
p.ylim((-9,0))
p.grid()
p.show()

p.savefig("M200c-cumulative-function-model.pdf")
p.clf()

#####################

fig = p.figure(1,(9,9))
ax = fig.add_subplot(111, projection='3d')
dd = ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5)
sc1 = ax.scatter(M200c, redshift,yData, s=n.ones_like(yData)*3, c='r', marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
ax.legend()
ax.set_xlabel(r'log $M_{200c}$ [km s$^{-1}$]')
ax.set_ylabel('redshift')
ax.set_ylim((0,n.max(redshift)))
ax.set_zlabel(r'log N($>M_{200c}$) [ h$^3$ Mpc$^{-3}$]')
ax.set_zlim((-8,0))
p.show()

#ax.set_yscale('log')
#ax.set_zscale('log')
p.savefig(join(Pdir , "M200c-cumulative-function-allZ-model.pdf"))
p.clf()

fig = p.figure(1,(9,9))
ax = fig.add_subplot(111, projection='3d')
sc1 = ax.scatter(M200c, redshift, 10**yData/10**vfGbis(M200c,redshift,res.x), s=n.ones_like(yData)*3, c='r', marker='o',label="MD data", rasterized=True)
sc1.set_edgecolor('face')
ax.legend()
ax.set_xlabel(r'log $M_{200c}$ [km s$^{-1}$]')
ax.set_ylabel('redshift')
ax.set_ylim((0,zmax))
ax.set_zlabel(r'Data / Model')
ax.set_zlim((0.5,1.5))
#ax.set_yscale('log')
#ax.set_zscale('log')
p.savefig(join(Pdir , "M200c-cumulative-function-allZ-modelRatio.pdf"))
p.clf()



sys.exit()
################################ Plot differential halo mass function  ################################


xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join(dir_40, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_4Gpc.dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join(dir_10, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_1Gpc"+".dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join(dir_25, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_2.5Gpc"+".dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join(dir_40, property_dir, type + "-"+ cos +"-"+ qty  +"_ALL_differential_MD_4Gpc.dat"),unpack=True)

#rhom_04 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value/aa.h**2 for zz in z_04])
#rhom_10 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value/aa.h**2 for zz in z_10])
#rhom_25 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value/aa.h**2 for zz in z_25])
#rhom_40 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value/aa.h**2 for zz in z_40])

rhom_04 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value for zz in z_04])
rhom_10 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value for zz in z_10])
rhom_25 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value for zz in z_25])
rhom_40 = n.array([aa.critical_density(zz).to(uu.solMass*uu.Mpc**-3).value for zz in z_40])

p.figure(0,(6,6))
p.axes([0.17,0.17,0.75,0.75])

#p.plot(n.log10(dat[0]),n.log10(dat[0]*dat[0]*dat[5]),'k--',lw=2)
s_04 = (z_04 <= 0.01)
y_04 = n.log10(yData_04*xData_04**2./rhom_04)
#y_04 = n.log10(yData_04*xData_04*xData_04)
p.plot(n.log10(xData_04[s_04]), y_04[s_04], marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)

s_10 = (z_10 == 0)
y_10 = n.log10(yData_10*xData_10**2./rhom_10)
#y_10 = n.log10(yData_10*xData_10*xData_10)
p.plot(n.log10(xData_10[s_10]),y_10[s_10], marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)

s_25 = (z_25 == 0)
y_25 = n.log10(yData_25*xData_25**2./rhom_25)
#y_25 = n.log10(yData_25*xData_25*xData_25)
p.plot(n.log10(xData_25[s_25]),y_25[s_25], marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)

s_40 = (z_40 == 0)
y_40 = n.log10(yData_40*xData_40**2./rhom_40)
#y_40 = n.log10(yData_40*xData_40*xData_40)
p.plot(n.log10(xData_40[s_40]),y_40[s_40], marker ='p', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)

p.xlabel(r'log$_{10}[M_{200c}/(h^{-1}M_\odot)]$')
p.ylabel(r'log$_{10}[(M^2/\rho_m) \; dn/dM]')
p.legend(loc=3)
p.grid()
p.savefig(join(Pdir , "M200c-diff-function-z0.pdf"))
p.clf()

fig = p.figure(1,(9,9))
ax = fig.add_subplot(111, projection='3d')

sc1 = ax.scatter(n.log10(xData_04),z_04,y_04, s=n.ones_like(z_04)*3, c='r', marker='o',label="SMD", rasterized=True)
sc1.set_edgecolor('face')

sc1 = ax.scatter(n.log10(xData_10),z_10,y_10, s=n.ones_like(z_10)*3, c='c', marker='v',label="MDPL", rasterized=True)
sc1.set_edgecolor('face')

sc1 = ax.scatter(n.log10(xData_25),z_25,y_25, s=n.ones_like(z_25)*3, c='m', marker='s',label="BigMD", rasterized=True)
sc1.set_edgecolor('face')

sc1 = ax.scatter(n.log10(xData_40),z_40,y_40, s=n.ones_like(z_40)*3, c='b', marker='p',label="HMD", rasterized=True)
sc1.set_edgecolor('face')

ax.legend()
ax.set_xlabel(r'$log_{10}[M_{200c}/(h^{-1}M_\odot)]$')
ax.set_ylabel('redshift')
ax.set_ylim((0,zmax))
ax.set_zlabel(r'$log_{10}[(M^2/\rho_m) \; dn/dM]$')
ax.set_zlim((-4,0))
p.savefig(join(Pdir , "M200c-diff-function-allZ.pdf"))
p.clf()

