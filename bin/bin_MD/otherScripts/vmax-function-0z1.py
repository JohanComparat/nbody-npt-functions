"""

ls Multidark-lightcones/MD_*/properties/vmax-mvir/hist-Central*1.0*.dat
ls Multidark-lightcones/MD_*/properties/vmax-mvir/hist-Central*0.9*.dat
ls Multidark-lightcones/MD_*/properties/vmax-mvir/hist-Central*0.8*.dat
ls Multidark-lightcones/MD_*/properties/vmax-mvir/hist-Central*0.7*.dat
ls Multidark-lightcones/MD_*/properties/vmax-mvir/hist-Central*0.6*.dat
ls Multidark-lightcones/MD_*/properties/vmax-mvir/hist-Central*0.5*.dat
ls Multidark-lightcones/MD_*/properties/vmax-mvir/hist-Central*0.4*.dat
ls Multidark-lightcones/MD_*/properties/vmax-mvir/hist-Central*0.3*.dat
ls Multidark-lightcones/MD_*/properties/vmax-mvir/hist-Central*0.2*.dat
"""
from scipy.interpolate import interp1d
import numpy as n
import matplotlib
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle
from os.path import join
from scipy.optimize import minimize

dir = "/Volumes/data/BigMD/vmaxFunction/"
Pdir = "/Volumes/data/BigMD/vmaxFunction/plots/"

def getVF(b0, b1, val,volume,label="SMD",completeness = 100, maxV=1500,errFactor=1.):
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = errFactor/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) &(yData_A * volume >= 1) &(xData_A > completeness)&(xData_A < maxV)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    return xData,yData,yDataErr,volume


def plotVFv3(b0, b1, val,volume,label="SMD",errFactor=1.):
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = errFactor/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) &(yData_A * volume >= 1)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    p.errorbar(xData,xData**3.1*yData,yerr=xData**3.1*yDataErr,elinewidth=1, label=label, rasterized=True)
    #print yData
    #p.axhline(1/(volume),label=r'1/V '+label,ls='dashed')


#saundersFct=lambda v, A, logv0,alpha,sig : 10**A * ( v /10**logv0)**(alpha) * n.e**( -n.log10( 1 + v /10**logv0)/(2*sig**2.))
#schechterFct=lambda v, A, logv0,alpha, sig : 10**A * ( v /10**logv0)**(alpha) * n.e**( - v / 10**logv0 /(2*sig**2.) )
#ps * (10**logl/10**logls)**(a+1) * n.e**(-10**logl/10**logls)
#doublePL=lambda v,A,logv0 ,a,b: 10**A * 10**( (1+a)*( v - 10**logv0) + b )



vf = lambda v, A, v0, alpha, beta : 10**A * (v/10**v0)**beta * n.e**(- (v/10**v0)**alpha )
vx = n.logspace(1.7,3.4,100)
pl = vf(vx,4,3,2.15,-2.9)

# for z<3
limits_04 = [80, 300]
limits_10 = [150, 400]
limits_25 = [400, 2000]
limits_40 = [800, 2500]
zmax = 1.

dir_04 = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/"
dir_10 = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/"
dir_25 = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/"
dir_40 = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/"

property_dir = "properties/vmax-mvir"
type = "hist"
cos = "Central"
qty = "vmax"
aa = "1.00000"
zz = 1./float(aa)-1.

fileName = type + "-"+ cos +"-"+ qty +"-"+aa+".dat"
fileName = type + "-"+ cos +"-"+ qty +"-*.dat"

fileList = glob.glob(join(dir_04, property_dir,fileName))
xData_04,yData_04,yDataErr_04,z_04 = [], [], [], []

for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
    xData_04_ii,yData_04_ii,yDataErr_04_ii,volume_04_ii = getVF(b0_04, b1_04, val_04,400.**3.,completeness = limits_04[0], maxV = limits_04[1])
    z_04_ii = (1/float(SMDfile.split('-')[-1][:-4])-1)*n.ones_like(xData_04_ii)
    if z_04_ii[0]<zmax :
        xData_04.append(xData_04_ii)
        yData_04.append(yData_04_ii)
        yDataErr_04.append(yDataErr_04_ii)
        z_04.append(z_04_ii)

z_04 = n.hstack((z_04))
xData_04 = n.hstack((xData_04))
yData_04 = n.hstack((yData_04))
yDataErr_04 = n.hstack((yDataErr_04))

n.savetxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_0.4Gpc"+".dat"),n.transpose([xData_04,z_04,yData_04,yDataErr_04]))

fileList = glob.glob(join(dir_10, property_dir,fileName))
xData_10,yData_10,yDataErr_10,z_10 = [], [], [], []

for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    b0_10, b1_10, val_10 = n.loadtxt(SMDfile,unpack=True)
    xData_10_ii,yData_10_ii,yDataErr_10_ii,volume_10_ii = getVF(b0_10, b1_10, val_10,1000.**3.,completeness = limits_10[0], maxV = limits_10[1])
    z_10_ii = (1/float(SMDfile.split('-')[-1][:-4])-1)*n.ones_like(xData_10_ii)
    if z_10_ii[0]<zmax :
        xData_10.append(xData_10_ii)
        yData_10.append(yData_10_ii)
        yDataErr_10.append(yDataErr_10_ii)
        z_10.append(z_10_ii)

z_10 = n.hstack((z_10))
xData_10 = n.hstack((xData_10))
yData_10 = n.hstack((yData_10))
yDataErr_10 = n.hstack((yDataErr_10))

n.savetxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_1Gpc"+".dat"),n.transpose([xData_10,z_10,yData_10,yDataErr_10]))


fileList = glob.glob(join(dir_25, property_dir,fileName))
xData_25,yData_25,yDataErr_25,z_25 = [], [], [], []

for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    b0_25, b1_25, val_25 = n.loadtxt(SMDfile,unpack=True)
    xData_25_ii,yData_25_ii,yDataErr_25_ii,volume_25_ii = getVF(b0_25, b1_25, val_25,2500.**3.,completeness = limits_25[0], maxV = limits_25[1])
    z_25_ii = (1/float(SMDfile.split('-')[-1][:-4])-1)*n.ones_like(xData_25_ii)
    if z_25_ii[0]<zmax :
        xData_25.append(xData_25_ii)
        yData_25.append(yData_25_ii)
        yDataErr_25.append(yDataErr_25_ii)
        z_25.append(z_25_ii)

z_25 = n.hstack((z_25))
xData_25 = n.hstack((xData_25))
yData_25 = n.hstack((yData_25))
yDataErr_25 = n.hstack((yDataErr_25))

n.savetxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_2.5Gpc"+".dat"),n.transpose([xData_25,z_25,yData_25,yDataErr_25]))


fileList = glob.glob(join(dir_40, property_dir,fileName))
xData_40,yData_40,yDataErr_40,z_40 = [], [], [], []

for ii in range(len(fileList)):
    SMDfile = fileList[ii]
    b0_40, b1_40, val_40 = n.loadtxt(SMDfile,unpack=True)
    xData_40_ii,yData_40_ii,yDataErr_40_ii,volume_40_ii = getVF(b0_40, b1_40, val_40,4000.**3.,completeness = limits_40[0], maxV = limits_40[1])
    z_40_ii = (1/float(SMDfile.split('-')[-1][:-4])-1)*n.ones_like(xData_40_ii)
    if z_40_ii[0]<zmax :
        xData_40.append(xData_40_ii)
        yData_40.append(yData_40_ii)
        yDataErr_40.append(yDataErr_40_ii)
        z_40.append(z_40_ii)

z_40 = n.hstack((z_40))
xData_40 = n.hstack((xData_40))
yData_40 = n.hstack((yData_40))
yDataErr_40 = n.hstack((yDataErr_40))

n.savetxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_4Gpc"+".dat"),n.transpose([xData_40,z_40,yData_40,yDataErr_40]))


xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_0.4Gpc"+".dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_1Gpc"+".dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_2.5Gpc"+".dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_4Gpc"+".dat"),unpack=True)

redshift = n.hstack(( z_04, z_10, z_25, z_40))
vmax = n.hstack(( xData_04, xData_10, xData_25, xData_40))
yData = n.hstack(( yData_04, yData_10, yData_25, yData_40))
yDataErr = n.hstack(( yDataErr_04, yDataErr_10, yDataErr_25, yDataErr_40))

vf = lambda v, A, v0, alpha, beta : 10**A * (v/10**v0)**beta * n.e**(- (v/10**v0)**alpha )

vfG = lambda v, z, A0, A1, vcut0, vcut1, a0, a1, a2, b0, b1 : 10**(A0 + A1 * z) * (v/10**(vcut0 + vcut1 * z))**(b0 + b1 * z) * n.e**(- (v/10**(vcut0 + vcut1 * z))**(a0 + a1 * z + a2 * z**2.) )

p1 = n.array([-3.48089246,  0.73575868,  2.73872233, -0.17861551,  1.47758206,  -0.32736035,  0.221566  , -2.70150905,  0.16587088])#4, 0.6, 2.85, -0.1, 1.77, -0.03, 0., -2.8, -0.1 ])

vfGbis = lambda v, z, p0 : vfG(v,z,p0[0],p0[1],p0[2],p0[3],p0[4],p0[5],p0[6],p0[7],p0[8])

chi2fun = lambda p0 : n.sum((vfGbis(vmax,redshift,p0) - yData)**2. / yDataErr**2. )/len(yDataErr)

print "looks for the optimum parameters"
res = minimize(chi2fun, p1, method='Powell',options={'xtol': 1e-6, 'disp': True, 'maxiter' : 50000000, 'nfev': 1800000})

print "ndof=",len(yDataErr)
print res
A0, A1, vcut0, vcut1, a0, a1, a2, b0, b1 = n.round(res.x,2)
str(A0)+"+"+str(A1)

print "A(z) & = "+str(A0)+" + "+str(A1)+r'\times z \\'
print r" V_{cut}(z) & = "+str(vcut0)+" + "+str(vcut1)+r'\times z \\'
print r" \alpha(z) & = "+str(a0)+" + "+str(a1)+r'\times z + '+str(a2)+r'\times z^2 \\'
print r" \beta(z) & = "+str(b0)+" + "+str(b1)+r'\times z \\'

# now outputs the model
X,Y = n.meshgrid(n.logspace(n.log10(limits_04[0]),n.log10(limits_40[1]),200),n.arange(0,zmax,0.025))

Z = vfGbis(X,Y,res.x)

n.savetxt(join("/Volumes/data/BigMD/vmaxFunction/data/best_fit.dat"),n.transpose([n.hstack((X)), n.hstack((Y)), n.hstack((Z))]) )

#######################################################
# now plots the results of the fit
print "now plots the results of the fit"

xData_04,z_04,yData_04,yDataErr_04 = n.loadtxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_0.4Gpc"+".dat"),unpack=True)
xData_10,z_10,yData_10,yDataErr_10 = n.loadtxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_1Gpc"+".dat"),unpack=True)
xData_25,z_25,yData_25,yDataErr_25 = n.loadtxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_2.5Gpc"+".dat"),unpack=True)
xData_40,z_40,yData_40,yDataErr_40 = n.loadtxt(join("/Volumes/data/BigMD/vmaxFunction/data/", type + "-"+ cos +"-"+ qty  +"MD_4Gpc"+".dat"),unpack=True)

vmax_mod, z_mod, n_mod = n.loadtxt(join("/Volumes/data/BigMD/vmaxFunction/data/best_fit.dat"), unpack=True)

#####################
from mpl_toolkits.mplot3d import Axes3D

fig = p.figure(1,(9,9))
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(n.log10(X), Y, n.log10(Z), rstride=10, cstride=10)

sc1 = ax.scatter(n.log10(xData_04),z_04,n.log10(yData_04), s=n.ones_like(z_04)*3, c='r', marker='o',label="SMD")
sc1.set_edgecolor('face')

sc1 = ax.scatter(n.log10(xData_10),z_10,n.log10(yData_10), s=n.ones_like(z_10)*3, c='c', marker='v',label="MDPL")
sc1.set_edgecolor('face')

sc1 = ax.scatter(n.log10(xData_25),z_25,n.log10(yData_25), s=n.ones_like(z_25)*3, c='m', marker='s',label="BigMD")
sc1.set_edgecolor('face')

sc1 = ax.scatter(n.log10(xData_40),z_40,n.log10(yData_40), s=n.ones_like(z_40)*3, c='b', marker='p',label="HMD")
sc1.set_edgecolor('face')

ax.legend()
ax.set_xlabel(r'log $V_{max}$ [km s$^{-1}$]')
ax.set_ylabel('redshift')
ax.set_ylim((0,1))
ax.set_zlabel(r'log N($>V_{max}$) [ h$^3$ Mpc$^{-3}$]')
#ax.set_yscale('log')
#ax.set_zscale('log')
p.savefig(Pdir + "vmax-function-3d-0-z-1.pdf")
p.show()



sys.exit()

######################
# z=0
zPlot=0.
s04 = (z_04==zPlot)
s10 = (z_10==zPlot)
s25 = (z_25==zPlot)
s40 = (z_40==zPlot)

p1 = n.array([-3.51218913,  0.77575315,  2.73736166, -0.17556627,  1.39837135,  -0.23156403,  0.262158  , -2.78723582,  0.31995556])
vfG = lambda v, z, A0, A1, vcut0, vcut1, a0, a1, a2, b0, b1 : 10**(A0 + A1 * z) * (v/10**(vcut0 + vcut1 * z))**(b0 + b1 * z) * n.e**(- (v/10**(vcut0 + vcut1 * z))**(a0 + a1 * z + a2 * z**2.) )

vfGbis = lambda v, z, p0 : vfG(v,z,p0[0],p0[1],p0[2],p0[3],p0[4],p0[5],p0[6],p0[7],p0[8])


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])

p.plot(xData_04[s04],yData_04[s04]*xData_04[s04]**2.8, marker ='o', mfc='None',mec='r',ls='none', label="SMD", rasterized=True)
p.plot(xData_10[s10],yData_10[s10]*xData_10[s10]**2.8, marker ='v', mfc='None',mec='c',ls='none', label="MDPL", rasterized=True)
p.plot(xData_25[s25],yData_25[s25]*xData_25[s25]**2.8, marker ='s', mfc='None',mec='m',ls='none', label="BigMD", rasterized=True)
p.plot(xData_40[s40],yData_40[s40]*xData_40[s40]**2.8, marker ='p', mfc='None',mec='b',ls='none', label="HMD", rasterized=True)

vmax_mod = n.hstack(( xData_04[s04], xData_10[s10], xData_25[s25], xData_40[s40] ))
vmax_mod.sort()
n_mod = vfGbis(vmax_mod, zPlot, p1)

p.plot(vmax_mod, n_mod* vmax_mod**2.8, 'k--', lw=2, label=r'best fit')

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) $V_{max}^{2.8}$ [ h$^3$ Mpc$^{-3}$ (km s$^{-1}$)$^{2.8}$]')
p.xscale('log')
p.xlim((80,3000))
p.yscale('log')
p.legend(loc=3)
p.grid()
p.savefig(Pdir + "vmax-function-z-"+str(zPlot)+".pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.9,1.1))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z0.00.pdf")
p.clf()

print n.round(z_25,3), " & ", n.round(res,3), "\\\\"



sys.exit()






























SMD_VFfile = join(dir_04, property_dir,fileName)
MDPL_VFfile = join(dir_10, property_dir,fileName)
BigMD_VFfile = join(dir_25, property_dir,fileName)
HMD_VFfile = join(dir_40, property_dir,fileName)


SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax-mvir/hist2d-Central-vmax-mvir-"+aa+".dat"
mvirBinsFile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax-mvir/mvir.bins"
vmaxBinsFile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax-mvir/vmax.bins"

MDPLfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax-mvir/hist2d-Central-vmax-mvir-"+aa+".dat"
BigMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/properties/vmax-mvir/hist2d-Central-vmax-mvir-"+aa+".dat"
HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax-mvir/hist2d-Central-vmax-mvir-"+aa+".dat"





H_04 = n.loadtxt(SMDfile)
mvir_bins_04 = n.loadtxt(mvirBinsFile)
vmax_bins_04 = n.loadtxt(vmaxBinsFile)
X,Y = n.meshgrid(( mvir_bins_04[1:] + mvir_bins_04[:-1])/2., (vmax_bins_04[1:]+vmax_bins_04[:-1])/2.)

sel = (H_04>0)

n.sum(H_04, axis=0)[n.sum(H_04, axis=0)>0]
n.sum(H_04, axis=1)[n.sum(H_04, axis=1)>0]

relation = interp1d(( mvir_bins_04[1:] + mvir_bins_04[:-1])/2.,H_04[180])
H_04[180].sum()



p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])
H_04 = n.loadtxt(SMDfile)
p.contour(X,Y,n.log10(H_04))#, colors='r')
H_10 = n.loadtxt(MDPLfile)
p.contour(X,Y,n.log10(H_10))#, colors='r')
#H_25 = n.loadtxt(BigMDfile)
#p.contour(X,Y,H_25)
H_40 = n.loadtxt(HMDfile)
p.contour(X,Y,n.log10(H_40))#, colors='r')
p.xlabel(r'log$_{10}(M_{vir})$')
p.ylabel(r'log$_{10}(V_{max})$')
p.colorbar()
p.savefig(Pdir + "mvir-vmax-z0.00.pdf")
p.show()


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 150)#, errFactor=10.)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

b0_10, b1_10, val_10 = n.loadtxt(MDPLfile,unpack=True)
z_10 = 1/float(MDPLfile.split('-')[-1][:-4])-1
#xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,label="MDPL",completeness = 0, maxV = 5000)
#p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,completeness = 151, maxV = 300)#, errFactor=10.)
p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=2,fmt='none', label="MDPL z="+str(n.round(z_10,3)), rasterized=True)

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 301, maxV = 500)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)

b0_40, b1_40, val_40 = n.loadtxt(HMDfile,unpack=True)
z_40 = 1/float(HMDfile.split('-')[-1][:-4])-1
#xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,label="HMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=1, ecolor='k',fmt='none', rasterized=True)
xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,completeness = 501, maxV = 2500)
p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=2,fmt='none', label="HMD z="+str(n.round(z_40,3)), rasterized=True)
p.axhline(1/(volume_40),label=r'1/V(HMD)',ls='dashed')

xData = n.hstack((xData_04, xData_10, xData_25, xData_40))
yData = n.hstack((yData_04, yData_10, yData_25, yData_40))
yDataErr = n.hstack((yDataErr_04, yDataErr_10, yDataErr_25, yDataErr_40))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)
#res, cov = curve_fit(saundersFct, xData, yData, sigma = yDataErr, p0 = (-3,2.5 ,-3.,1.2) , maxfev = 5000000)
#res, cov = curve_fit(schechterFct, xData, yData, sigma = yDataErr, p0 = (-3,2.5 ,-3.,1.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
#p.plot(vx, saundersFct(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
#p.plot(vx, schechterFct(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.xlim((40,3000))
p.ylim((0.5/(volume_40),1))
p.yscale('log')
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z0.00.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.9,1.1))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z0.00.pdf")
p.clf()

print n.round(z_25,3), " & ", n.round(res,3), "\\\\"



######### z = 0.05 a = 0.950 ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.95600.dat"
MDPLfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.95670.dat"
BigMDfile = "/Volumes/data/BigMD/2.5Gpc/vmax/hist-Central-vmax-0.95600.dat"
#HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax/hist-Central-vmax-0.25320.dat"


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 220, maxV = 2500)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)
p.axhline(1/(volume_25),label=r'1/V(BigMD)',ls='dashed')

b0_10, b1_10, val_10 = n.loadtxt(MDPLfile,unpack=True)
z_10 = 1/float(MDPLfile.split('-')[-1][:-4])-1
#xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,label="MDPL",completeness = 0, maxV = 5000)
#p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,completeness = 120, maxV = 220)
p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=2,fmt='none', label="MDPL z="+str(n.round(z_10,3)), rasterized=True)

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 120)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

xData = n.hstack((xData_04, xData_10, xData_25))
yData = n.hstack((yData_04, yData_10, yData_25))
yDataErr = n.hstack((yDataErr_04, yDataErr_10, yDataErr_25))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.ylim((0.5/(volume_40),1))
p.xlim((40,3000))
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z0.05.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.95,1.05))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z0.05.pdf")
p.clf()

print n.round(z_25,3), " & ", n.round(res,3), "\\\\"



######### z = 0.42 a = 0.70 ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.70030.dat"
MDPLfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.70160.dat"
BigMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/properties/vmax/hist-Central-vmax-0.70030.dat"
#HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax/hist-Central-vmax-0.25320.dat"


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 220, maxV = 2500)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)
p.axhline(1/(volume_25),label=r'1/V(BigMD)',ls='dashed')

b0_10, b1_10, val_10 = n.loadtxt(MDPLfile,unpack=True)
z_10 = 1/float(MDPLfile.split('-')[-1][:-4])-1
#xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,label="MDPL",completeness = 0, maxV = 5000)
#p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,completeness = 120, maxV = 220)
p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=2,fmt='none', label="MDPL z="+str(n.round(z_10,3)), rasterized=True)

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 120)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

xData = n.hstack((xData_04, xData_10, xData_25))
yData = n.hstack((yData_04, yData_10, yData_25))
yDataErr = n.hstack((yDataErr_04, yDataErr_10, yDataErr_25))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.ylim((0.5/(volume_40),1))
p.xlim((40,3000))
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z0.42.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.95,1.05))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z0.42.pdf")
p.clf()

print n.round(z_25,3), " & ", n.round(res,3), "\\\\"


######### z = 0.88 a = 0.5 ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.53000.dat"
BigMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/properties/vmax/hist-Central-vmax-0.53000.dat"
HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax/hist-Central-vmax-0.53780.dat"


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])
b0_40, b1_40, val_40 = n.loadtxt(HMDfile,unpack=True)
z_40 = 1/float(HMDfile.split('-')[-1][:-4])-1
#xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,label="HMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=1, ecolor='k',fmt='none', rasterized=True)
xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,completeness = 1000, maxV = 2500)
p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=2,fmt='none', label="HMD z="+str(n.round(z_40,3)), rasterized=True)
p.axhline(1/(volume_40),label=r'1/V(HMD)',ls='dashed')

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 200, maxV = 1000)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 200)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

xData = n.hstack((xData_04, xData_25, xData_40))
yData = n.hstack((yData_04, yData_25, yData_40))
yDataErr = n.hstack((yDataErr_04, yDataErr_25, yDataErr_40))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.xlim((40,3000))
p.ylim((0.5/(volume_40),1))
p.yscale('log')
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z0.88.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.95,1.05))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z0.88.pdf")
p.clf()


print n.round(z_25,3), " & ", n.round(res,3), "\\\\"



######### z = 3 a = 0.25 ############

SMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.24800.dat"
MDPLfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.24770.dat"
BigMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_2.5Gpc/properties/vmax/hist-Central-vmax-0.25700.dat"
HMDfile = "/Volumes/data/BigMD/Multidark-lightcones/MD_4Gpc/properties/vmax/hist-Central-vmax-0.25320.dat"


p.figure(1,(6,6))
p.axes([0.17,0.17,0.75,0.75])
b0_40, b1_40, val_40 = n.loadtxt(HMDfile,unpack=True)
z_40 = 1/float(HMDfile.split('-')[-1][:-4])-1
#xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,label="HMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=1, ecolor='k',fmt='none', rasterized=True)
xData_40,yData_40,yDataErr_40,volume_40 = getVF(b0_40, b1_40, val_40,4000.**3.,completeness = 1030, maxV = 2500)
p.errorbar(xData_40,yData_40,yerr=yDataErr_40,elinewidth=2,fmt='none', label="HMD z="+str(n.round(z_40,3)), rasterized=True)
p.axhline(1/(volume_40),label=r'1/V(HMD)',ls='dashed')

b0_25, b1_25, val_25 = n.loadtxt(BigMDfile,unpack=True)
z_25 = 1/float(BigMDfile.split('-')[-1][:-4])-1
#xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_25,yData_25,yDataErr_25,volume_25 = getVF(b0_25, b1_25, val_25, 2500.**3.,completeness = 220, maxV = 1030)
p.errorbar(xData_25,yData_25,yerr=yDataErr_25,elinewidth=2,fmt='none', label="BigMD z="+str(n.round(z_25,3)), rasterized=True)

b0_10, b1_10, val_10 = n.loadtxt(MDPLfile,unpack=True)
z_10 = 1/float(MDPLfile.split('-')[-1][:-4])-1
#xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,label="MDPL",completeness = 0, maxV = 5000)
#p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_10,yData_10,yDataErr_10,volume_10 = getVF(b0_10, b1_10, val_10, 1000.**3.,completeness = 120, maxV = 220)
p.errorbar(xData_10,yData_10,yerr=yDataErr_10,elinewidth=2,fmt='none', label="MDPL z="+str(n.round(z_10,3)), rasterized=True)

b0_04, b1_04, val_04 = n.loadtxt(SMDfile,unpack=True)
z_04 = 1/float(SMDfile.split('-')[-1][:-4])-1
#xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,label="SMD",completeness = 0, maxV = 5000)
#p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=1,ecolor='k',fmt='none', rasterized=True)
xData_04,yData_04,yDataErr_04,volume_04 = getVF(b0_04, b1_04, val_04,400.**3.,completeness = 50, maxV = 120)
p.errorbar(xData_04,yData_04,yerr=yDataErr_04,elinewidth=2,fmt='none', label="SMD z="+str(n.round(z_04,3)), rasterized=True)

xData = n.hstack((xData_04, xData_10, xData_25, xData_40))
yData = n.hstack((yData_04, yData_10, yData_25, yData_40))
yDataErr = n.hstack((yDataErr_04, yDataErr_10, yDataErr_25, yDataErr_40))
res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)

p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))

p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.ylim((0.5/(volume_40),1))
p.xlim((40,3000))
p.legend(loc=3)
p.grid()
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.savefig(Pdir + "fit-perBox-z3.0.pdf")
p.clf()


p.figure(2,(6,6))
p.axes([0.17,0.15,0.75,0.75])
p.fill_between(xData,1+yDataErr / yData,1-yDataErr / yData,rasterized=True,alpha=0.5)
p.plot(xData, yData / vf(xData,res[0], res[1], res[2], res[3]),'k')
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$) / best fit')
p.xscale('log')
p.ylim((0.9,1.1))
p.xlim((40,3000))
p.title("A="+str(n.round(res[0],2))+r" v$_0$="+ str(n.round(res[1],2))+r" $\alpha$="+ str(n.round(res[2],2))+r" $\beta$="+ str(n.round(res[3],2)) )
p.grid()
p.savefig(Pdir + "fit-ratio-z3.0.pdf")
p.clf()


print n.round(z_25,3), " & ", n.round(res,3), "\\\\"



















p.figure(1,(6,6))
p.axes([0.2,0.2,0.75,0.75])

b0_10, b1_10, val_10 = n.loadtxt("Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.24770.dat",unpack=True)
xData,yData,yDataErr,volume = plotVF(b0_10, b1_10, val_10,1000.**3.,label="MDPL",completeness = 100, maxV = 1500)

p.errorbar(xData,yData,yerr=yDataErr,elinewidth=1, label=label, rasterized=True)
p.axhline(1/(volume),label=r'1/V '+label,ls='dashed')

b0_25, b1_25, val_25 = n.loadtxt("Multidark-lightcones/MD_2.5Gpc/properties/vmax/hist-Central-vmax-0.25700.dat",unpack=True)
plotVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD")


plotVF(b0_25, b1_25, val_25, 4000.**3.,label="HMD")

p.xlim((50,4000))
p.ylim((0.5/(4000.**3.), 1))
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.legend(loc=3)
p.savefig(Pdir + "VmaxF-cumulative-central-z3.pdf")


p.figure(2,(6,6))
p.axes([0.2,0.2,0.75,0.75])

b0_04, b1_04, val_04 = n.loadtxt("Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.71830.dat",unpack=True)
plotVF(b0_04, b1_04, val_04,400.**3.,label="SMD")

b0_10, b1_10, val_10 = n.loadtxt("Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.71730.dat",unpack=True)
plotVF(b0_10, b1_10, val_10,1000.**3.,label="MDPL")

b0_25, b1_25, val_25 = n.loadtxt("2.5Gpc/vmax/hist-Central-vmax-0.71830.dat",unpack=True)
plotVF(b0_25, b1_25, val_25, 2500.**3.,label="BigMD")


p.xlim((5,4000))
p.ylim((0.5/(2500.**3.), 10))
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.legend(loc=3)
p.grid()
p.savefig(Pdir + "VmaxF-cumulative-central-z0.42.pdf")
p.clf()



p.figure(2,(6,6))
p.axes([0.2,0.2,0.75,0.75])

#b0_04, b1_04, val_04 = n.loadtxt("Multidark-lightcones/MD_0.4Gpc/properties/vmax/hist-Central-vmax-0.71830.dat",unpack=True)
#plotVFv3(b0_04, b1_04, val_04,400.**3.,label="SMD")

b0_10, b1_10, val_10 = n.loadtxt("Multidark-lightcones/MD_1Gpc_new_rockS/properties/vmax/hist-Central-vmax-0.71730.dat",unpack=True)
plotVFv3(b0_10, b1_10, val_10,1000.**3.,label="MDPL")

b0_25, b1_25, val_25 = n.loadtxt("2.5Gpc/vmax/hist-Central-vmax-0.71830.dat",unpack=True)
plotVFv3(b0_25, b1_25, val_25, 2500.**3.,label="BigMD")


p.xlim((100,400))
p.ylim((40000,50000))
p.xlabel(r'$V_{max}$ [km s$^{-1}$]')
p.ylabel(r'$V_{max}^{3.1}$N($>V_{max}$)  [ h$^3$ Mpc$^{-3}$ ]')
#p.xscale('log')
#p.yscale('log')
p.legend(loc=3)
p.grid()
p.savefig(Pdir + "VmaxFv3-cumulative-central-z0.42-zoom-o2.pdf")
p.clf()





centralList = n.array(glob.glob(dir+"hist-central-Vpeak-?.?????.dat"))
centralList.sort()

volume = (2500.)**3. 

# evolution plot

ids = [2,3,4,5,7,9,11,20,56]
p.figure(1,(6,6))
p.axes([0.2,0.2,0.75,0.75])
for iii in ids :
    el = centralList[iii]
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = 1/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) & (xData_A>20)&(yData_A * volume >= 1)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    p.errorbar(xData,yData,yerr=yDataErr,elinewidth=1, label="z="+str(n.round(z,3)), rasterized=True)

p.axhline(1/(volume),label=r'1/V',color='k',ls='dashed')
p.xlim((50,4000))
p.ylim((0.5/(volume), 1e-2))
p.xlabel(r'$V_{peak}$ [km s$^{-1}$]')
p.ylabel(r'N($>V_{peak}$)  [ h$^3$ Mpc$^{-3}$ ]')
p.xscale('log')
p.yscale('log')
p.legend(loc=3)
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z8-z0-evolution.pdf")
p.clf()

sys.exit()


filename = n.array(glob.glob("0.4Gpc/*/hist*Central-Vpeak-0.22*.dat"))[0]
def getData(filename,Vcut=0,volume = 400.**3.):
    zz=1./float(filename.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(filename,unpack=True)
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = 1/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) & (xData_A>Vcut)&(yData_A * volume >= 1)
    velocity = xData_A[sel]
    redshift = n.ones_like(velocity)*zz
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    return redshift, velocity, yData,yDataErr

z,v,y,yerr = getData(hist04)


hist10 = n.array(glob.glob("1Gpc/*/hist*Central-Vpeak-0.22*.dat"))[0]
z10=1./float(hist10.split('-')[-1][:-4])-1.
hist25 = n.array(glob.glob("2.5Gpc/*/hist*central-Vpeak-0.22*.dat"))[0]
z25=1./float(hist25.split('-')[-1][:-4])-1.



# cumulative velocity function

results, errors, redshifts, chi2Rs = [], [], [], []

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    Nc = n.array([ n.sum(val[ii:]) for ii in range(len(val)) ])
    xData_A = (10**b0+10**b1)/2.
    yData_A = Nc/(volume ) 
    yDataErr_A = 1/val**0.5
    sel = (yDataErr_A!=n.inf)&(yData_A>0) & (xData_A>200)&(yData_A * volume >= 1)
    xData = xData_A[sel]
    yData = yData_A[sel]
    yDataErr = yDataErr_A[sel]*yData_A[sel]
    res, cov = curve_fit(vf, xData, yData, sigma = yDataErr, p0 = (4,3,2.15,-3.2) , maxfev = 5000000)
    chi2red = n.round(n.sum( (vf(xData,res[0], res[1], res[2], res[3])- yData)**2. / (yDataErr**2)) / (len(xData) - (4-1)),2)
    chi2Rs.append(chi2red)
    results.append(res)
    errors.append(cov)
    redshifts.append(z)
    p.figure(1,(6,6))
    p.axes([0.2,0.2,0.75,0.75])
    p.plot(vx, vf(vx,res[0], res[1], res[2], res[3]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
    p.errorbar(xData,yData,yerr=yDataErr,color='k',elinewidth=2, label="z="+str(n.round(z,3)), rasterized=True)
    p.axhline(1/(volume),label=r'1/V',color='k',ls='dashed')
    p.xlim((100,4000))
    p.ylim((0.5/(volume), 1e-1))
    p.xlabel(r'$V_{peak}$ [km s$^{-1}$]')
    p.ylabel(r'N($>V_{peak}$)  [ h$^3$ Mpc$^{-3}$ ]')
    p.xscale('log')    
    p.yscale('log')    
    p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2))+" "+ str(n.round(res[2],2))+" "+ str(n.round(res[3],2)) )
    p.legend()
    p.grid()
    p.savefig(Pdir + "VpeakF-cumulative-central-z-"+str(n.round(z,4))+".pdf")
    p.clf()

results = n.transpose(results)
errors = n.array(errors)
diag_err = n.array([ [el[0][0]**0.5, el[1][1]**0.5, el[2][2]**0.5, el[3][3]**0.5] for el in errors]).T
redshifts = n.array(redshifts)
chi2Rs = n.array(chi2Rs)

f=open(Pdir + "VpeakF-cumulative-central-z-param-fit-results.pkl",'w')
cPickle.dump([results, errors, diag_err, redshifts, chi2Rs],f)
f.close()

ok = (redshifts < 1.1)

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[0][ok], diag_err[0][ok])
p.ylabel('log Amplitude')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-amplitude-evolution.pdf")
p.clf()
p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[1][ok], diag_err[1][ok])
p.ylabel(r'$V_0$')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-v0-evolution.pdf")
p.clf()
p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[2][ok], diag_err[2][ok])
p.ylabel(r'$\alpha$')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-alpha-evolution.pdf")
p.clf()
p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[3][ok], diag_err[3][ok])
p.ylabel(r'$\beta$')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-beta-evolution.pdf")
p.clf()
sys.exit()

Pdir = "/Volumes/data/BigMD/2.5Gpc/plots/"
dir = "/Volumes/data/BigMD/2.5Gpc/"

ARbins = n.loadtxt('/Volumes/data/BigMD/2.5Gpc/AccRate.bins')
Vbins = n.loadtxt('/Volumes/data/BigMD/2.5Gpc/Vpeak.bins')
Mbins = n.loadtxt('/Volumes/data/BigMD/2.5Gpc/Mpeak.bins')


centralList = n.array(glob.glob(dir+"hist2d-central-Mpeak-Accrate-?.?????.dat"))
centralList.sort()
print centralList

X,Y= n.meshgrid(ARbins, Mbins)

volume = (2500.)**3. 
norm = volume * n.median(Mbins[1:]-Mbins[:-1]) *n.median(ARbins[1:]-ARbins[:-1]) 
print norm , n.log10(0.8/norm)

p.figure(1,(10,6))
for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    data = n.loadtxt(el,unpack=True)
    p.title("z="+str(n.round(z,3)))
    p.pcolormesh(X, Y, n.log10(data.T/norm),vmin=n.log10(0.8/norm),vmax=-4,rasterized=True)
    cb = p.colorbar(shrink=0.7)
    cb.set_label('N/(V dlogM dAR) ')
    p.ylim((11,16))
    p.ylabel(r'log $M_{peak}$ [km/s]')
    p.xlim((0,50000))
    p.xlabel('Accretion rate Msun/h/yr at Mpeak')
    p.savefig(Pdir + "plot-hist2d-central-Mpeak-Accrate-"+str(n.round(z,5))+".pdf")
    p.clf()


sys.exit()



centralList = n.array(glob.glob(dir+"hist2d-central-Vpeak-Accrate-?.?????.dat"))
centralList.sort()
print centralList

X,Y= n.meshgrid(ARbins, Vbins)

volume = (2500.)**3. 
norm = volume * n.median(Vbins[1:]-Vbins[:-1]) *n.median(ARbins[1:]-ARbins[:-1]) 
p.figure(1,(10,6))
for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    data = n.loadtxt(el,unpack=True)
    p.title("z="+str(n.round(z,3)))
    p.pcolormesh(X, Y, n.log10(data.T/norm),vmin=n.log10(0.8/norm),vmax=-4,rasterized=True)
    cb = p.colorbar(shrink=0.7)
    cb.set_label('N/(V dlogV dAR) ')
    p.ylim((1.5,3.5))
    p.ylabel(r'log $V_{peak}$ [km/s]')
    p.xlim((0,50000))
    p.xlabel('Growth Rate of Mpeak [z, z+0.5] Msun/h/yr')
    p.savefig(Pdir + "plot-hist2d-central-Vpeak-Accrate-"+str(n.round(z,5))+".pdf")
    p.clf()





centralList = n.array(glob.glob(dir+"hist-sat-Mpeak-?.?????.dat"))
centralList.sort()
print centralList

p.figure(1,(10,6))

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    p.plot((10**b0+10**b1)/2.,val/(volume * (b1 - b0) ),label="z="+str(n.round(z,3)), rasterized=True)
    
    if iii%15==0 :
        p.axhline(1/(volume*n.median((b1 - b0))),label=r'1/(dlogMx2500$^3$)')
        p.axvline(23593750000.0*100, label='100 Mp')
        p.xlim((23593750000.0*50,5e16))
        p.ylim((0.5/(volume*n.median((b1 - b0))), 1e-1))
        p.xlabel(r'$M_{peak}$ [M$_\odot$/h]')
        p.ylabel(r'N($V_1<M_{peak}<V_2$) / dlog Mpeak / Volume [ h/M$_\odot$ . h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-sat-Mpeak-"+str(iii)+".pdf")
        p.clf()



centralList = n.array(glob.glob(dir+"hist-central-Mpeak-?.?????.dat"))
centralList.sort()
print centralList

volume = (2500.)**3. 

p.figure(1,(10,6))

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    p.plot((10**b0+10**b1)/2.,val/(volume * (b1 - b0) ),label="z="+str(n.round(z,3)), rasterized=True)

    if iii%15==0 :
        p.axhline(1/(volume*n.median((b1 - b0))),label=r'1/(dlogMx2500$^3$)')
        p.axvline(23593750000.0*100, label='100 Mp')
        p.xlim((23593750000.0*50,5e16))
        p.ylim((0.5/(volume*n.median((b1 - b0))), 1e-1))
        p.xlabel(r'$M_{peak}$ [M$_\odot$/h]')
        p.ylabel(r'N($V_1<M_{peak}<V_2$) / dlog Mpeak / Volume [ h/M$_\odot$. h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-central-Mpeak-"+str(iii)+".pdf")
        p.clf()


sys.exit()





volume = (2500.)**3. 

centralList = n.array(glob.glob(dir+"hist-sat-Vpeak-?.?????.dat"))
centralList.sort()


p.figure(1,(10,6))

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    p.plot((10**b0+10**b1)/2.,val/(volume * (b1 - b0) ),label="z="+str(n.round(z,3)), rasterized=True)
    
    if iii%15==0 :
        p.axhline(1/(volume*n.median((b1 - b0))),label=r'1/(0.01x2500$^3$)')
        p.xlim((100,1e5))
        p.ylim((0.5/(volume*n.median((b1 - b0))), 1e-1))
        p.xlabel(r'$V_{peak}$ [km/s]')
        p.ylabel(r'N($V_1<V_{peak}<V_2$) / dlog Vpeak / Volume [ s/km . h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-sat-Vpeak-"+str(iii)+".pdf")
        p.clf()



centralList = n.array(glob.glob(dir+"hist-central-Vpeak-?.?????.dat"))
centralList.sort()

volume = (2500.)**3. 

p.figure(1,(10,6))

for iii,el in enumerate(centralList) :
    z=1./float(el.split('-')[-1][:-4])-1.
    b0,b1,val = n.loadtxt(el,unpack=True)
    p.plot((10**b0+10**b1)/2.,val/(volume * (b1 - b0) ),label="z="+str(n.round(z,3)), rasterized=True)

    if iii%15==0 :
        p.axhline(1/(volume*n.median((b1 - b0))),label=r'1/(0.01x2500$^3$)')
        p.xlim((100,1e5))
        p.ylim((0.5/(volume*n.median((b1 - b0))), 1e-1))
        p.xlabel(r'$V_{peak}$ [km/s]')
        p.ylabel(r'N($V_1<V_{peak}<V_2$) / dlog Vpeak / Volume [ s/km . h3/Mpc3 ]')
        p.xscale('log')    
        p.yscale('log')    
        p.legend()
        p.grid()
        p.savefig(Pdir + "plot-hist-central-Vpeak-"+str(iii)+".pdf")
        p.clf()


sys.exit()




zzz=[]
for ii in range(len(centralList)):
    zzz.append( centralList[ii].split('-')[-1][:-4] )

Zset = list(set(zzz))

for el in Zset :
    centralList = n.array(glob.glob(dir+"hist-central-Vpeak-"+el+".dat"))
    centralList.sort()


    massB[:-1], massB[1:],  nnM.sum(axis=0)


    nnM = n.empty( [len(centralList),len(massB)-1] ) 
    nnV = n.empty( [len(centralList),len(vcirB)-1] )
    dataVC = n.empty( [len(centralList),len(vcirB)-1,len(concB)-1] )
    dataMC = n.empty( [len(centralList),len(massB)-1,len(concB)-1] )

    for jj in range(len(centralList)):
        f=open(centralList[jj],'r')
        nnMinter,nnVinter,nnCinter,dataMCinter,dataVCinter = cPickle.load(f)
        nnM[jj] = nnMinter
        nnV[jj] = nnVinter 
        dataMC[jj] = dataMCinter[0]
        dataVC[jj] = dataVCinter[0]
        f.close()


    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-central-Mpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([massB[:-1], massB[1:],  nnM.sum(axis=0)]))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-central-Vpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([vcirB[:-1], vcirB[1:],  nnV.sum(axis=0)]) )

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-central-Mpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataMC.sum(axis=0))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-central-Vpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataVC.sum(axis=0))


    satList = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/*"+box.get_snl()[ii].split('_')[-1][:-5]+"*MVAmatrixSatellite*.pkl"))
    satList.sort()

    nnM = n.empty( [len(satList),len(massB)-1] ) 
    nnV = n.empty( [len(satList),len(vcirB)-1] )
    dataVC = n.empty( [len(satList),len(vcirB)-1,len(concB)-1] )
    dataMC = n.empty( [len(satList),len(massB)-1,len(concB)-1] )

    for jj in range(len(satList)):
        f=open(satList[jj],'r')
        nnMinter,nnVinter,nnCinter,dataMCinter,dataVCinter = cPickle.load(f)
        nnM[jj] = nnMinter
        nnV[jj] = nnVinter 
        dataMC[jj] = dataMCinter[0]
        dataVC[jj] = dataVCinter[0]
        f.close()


    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-sat-Mpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([massB[:-1], massB[1:],  nnM.sum(axis=0)]))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-sat-Vpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([vcirB[:-1], vcirB[1:],  nnV.sum(axis=0)]) )

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-sat-Mpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataMC.sum(axis=0))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-sat-Vpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataVC.sum(axis=0))



