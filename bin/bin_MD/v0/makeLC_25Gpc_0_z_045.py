from lightcone import *


wdir = "/home2/jcomparat/eBOSS-LC/Multidark-lightcones/"
boxDir = "MD_2.5Gpc/"
snlInter=n.array(glob.glob(wdir+boxDir+"h*.list"))
print snlInter
zInter=n.array([ 1/float(el.split('_')[-1][:-5])-1 for el in snlInter ])
ids=n.argsort(zInter)
snl=snlInter[ids]
print snl
zs=zInter[ids]
zArray=n.arange(0.,0.5,1e-4)
Hbox=67.77 * uu.km / (uu.s * uu.Mpc)
outputName = wdir + boxDir + "MD_2.5Gpc_lc_square_0z0.5.txt"
Melement = 23593750000.0

lc1=LightconeSquare(zMean=0,Lbox=2500.0 * uu.Mpc,wdir=wdir,boxDir=boxDir,snl=snl,zsl=zs,zArray=zArray,Hbox=Hbox,outputName=outputName,Melement = Melement)
print "lc object created, defines limits now: "
lc1.defineLClimitsObsInTheCenter(massFactor=50)
print "limits ok"
lc1.defineSnapshotLimits()
print "starts lc building"
lc1.constructLC()




