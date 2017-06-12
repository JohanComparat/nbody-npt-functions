
import numpy as n
import matplotlib
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as p
import glob
import sys
from scipy.optimize import curve_fit
import cPickle

Pdir = "/Volumes/data/BigMD/2.5Gpc/propertiesAtPeak/plots/"
f=open(Pdir + "VpeakF-cumulative-central-z-param-fit-results.pkl",'r')
results, errors, diag_err, redshifts, chi2Rs = cPickle.load(f)
f.close()

Ns=3
"""
id = n.argsort(redshifts)
for ii in id : #range(id): 2))#len(redshifts)):
    print n.round(redshifts[ii],Ns)," & ", n.round(results[0][ii],Ns), r"$\pm$ ", n.round(diag_err[0][ii],Ns)," & ", n.round(results[1][ii],Ns),r" $\pm$ ", n.round(diag_err[1][ii],Ns)," & ", n.round(results[2][ii],Ns),r"$\pm$ ", n.round(diag_err[2][ii],Ns)," & ", n.round(results[3][ii],Ns)," $\pm$ ", n.round(diag_err[3][ii],Ns)," \\\\"#, n.round(chi2Rs[ii],Ns)," \\\\"

"""

ok = (redshifts < 1.1)

pl1 = lambda x, a, b : a * x + b
pl2 = lambda x, a0, a1, a2 : a0+ a1 * x + a2*x*x

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[0][ok], diag_err[0][ok])
res, cov = curve_fit(pl1, redshifts[ok],results[0][ok], sigma =  diag_err[0][ok], p0 = (-0.1,5) , maxfev = 5000000)
p.plot(redshifts[ok], pl1(redshifts[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.ylabel('log Amplitude')
p.xlabel('z')
p.grid()
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))
p.legend()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-amplitude-evolution.pdf")
p.clf()
print "A=(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),")"

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[1][ok], diag_err[1][ok])
res, cov = curve_fit(pl1, redshifts[ok],results[1][ok], sigma =  diag_err[1][ok], p0 = (-0.5,3) , maxfev = 5000000)
p.plot(redshifts[ok], pl1(redshifts[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))
p.legend()
p.ylabel(r'$V_0$')
p.xlabel('z')
p.grid()
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-v0-evolution.pdf")
p.clf()
print "V_0=(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),")"



p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[2][ok], diag_err[2][ok])
res, cov = curve_fit(pl2, redshifts[ok],results[2][ok], sigma =  diag_err[2][ok], p0 = (-0.7,2,1) , maxfev = 5000000)
p.plot(redshifts[ok], pl2(redshifts[ok],res[0], res[1],res[2]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.legend()
p.ylabel(r'$\alpha$')
p.xlabel('z')
p.grid()
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2))+" "+ str(n.round(res[2],2)))#+" "+ str(n.round(res[3],2)) )
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-alpha-evolution.pdf")
p.clf()
print "\\alpha=(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") z + (", n.round(res[2],Ns), r"\pm ", n.round(cov[2][2]**0.5,Ns),")z^2"

p.figure(0,(6,6))
p.axes([0.2,0.2,0.75,0.75])
p.errorbar(redshifts[ok],results[3][ok], diag_err[3][ok])
res, cov = curve_fit(pl1, redshifts[ok],results[3][ok], sigma =  diag_err[3][ok], p0 = (-0.5,3) , maxfev = 5000000)
p.plot(redshifts[ok], pl1(redshifts[ok],res[0], res[1]), label=r'fit')#$\chi^2$/ndof='+str(chi2red))
p.legend()
p.ylabel(r'$\beta$')
p.xlabel('z')
p.grid()
p.title(str(n.round(res[0],2))+" "+ str(n.round(res[1],2)))#+" "+ str(n.round(res[2],2))+" "+ str(n.round(res[3],2)) )
p.savefig(Pdir + "VpeakF-cumulative-central-z-param-beta-evolution.pdf")
p.clf()
print "\\beta=(", n.round(res[0],Ns), r"\pm ", n.round(cov[0][0]**0.5,Ns), ") z + (", n.round(res[1],Ns), r"\pm ", n.round(cov[1][1]**0.5,Ns),") "

