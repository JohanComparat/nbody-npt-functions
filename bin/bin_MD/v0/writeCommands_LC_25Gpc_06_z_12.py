# write all the python scripts to be submitted

from lightconeParallel import *

wdir = "/home2/jcomparat/eBOSS-LC/Multidark-lightcones/"
boxDir = wdir + "MD_2.5Gpc/"
lcDir = boxDir + "lc_square_0.6z1.2/"
snlInter=n.array(glob.glob(boxDir+"h*.list"))
print snlInter
zInter=n.array([ 1/float(el.split('_')[-1][:-5])-1 for el in snlInter ])
ids=n.argsort(zInter)
snl=snlInter[ids]
print snl
zs=zInter[ids]
zArray=n.arange(0.1,2.4,1e-4)
Hbox=67.77 * uu.km / (uu.s * uu.Mpc)
outputName = lcDir + "snapshot.txt"
Melement = 23593750000.0

lc1=LightconeSquare(zMean=0.88,Lbox=2500.0 * uu.Mpc,wdir=wdir,boxDir=boxDir,lcDir=lcDir ,snl=snl, zsl=zs,zArray=zArray,Hbox=Hbox,outputName=outputName,Melement = Melement)
print "lc object created, defines limits now: "
lc1.defineLClimits(massFactor=50)
print "limits ok"
lc1.defineSnapshotLimits()
print "starts lc building"


for ii in range(len(lc1.get_snl())):
	name=lc1.get_snl()[ii].split("/")[-1]
	f=open("run_"+name+'.py','a')
	f.write("""from lightconeParallel import *
wdir = "/home2/jcomparat/eBOSS-LC/Multidark-lightcones/"
boxDir = wdir + "MD_2.5Gpc/"
lcDir = boxDir + "lc_square_0.6z1.2/"
snlInter=n.array(glob.glob(boxDir+"h*.list"))
print snlInter
zInter=n.array([ 1/float(el.split('_')[-1][:-5])-1 for el in snlInter ])
ids=n.argsort(zInter)
snl=snlInter[ids]
print snl
zs=zInter[ids]
zArray=n.arange(0.1,2.4,1e-4)
Hbox=67.77 * uu.km / (uu.s * uu.Mpc)
outputName = lcDir + "snapshot.txt"
Melement = 23593750000.0
lc1=LightconeSquare(zMean=0.88,Lbox=2500.0 * uu.Mpc,wdir=wdir,boxDir=boxDir,lcDir=lcDir,snl=snl,zsl=zs,zArray=zArray,Hbox=Hbox,outputName=outputName,Melement = Melement)
print "lc object created, defines limits now: "
lc1.defineLClimits(massFactor=50)
print "limits ok"
lc1.defineSnapshotLimits()
print "starts lc building" \n""")
	f.write("""lc1.set_outputName(outputName = lcDir + "shell_" + '"""+ name +"""') \n""")
	f.write("""lc1.constructLC_snapshot("""+str(ii)+""")""")
	f.close()



