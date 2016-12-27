import fileinput
import astropy.cosmology as co
c2=co.Planck13
from scipy.interpolate import interp1d

import astropy.units as uu
import numpy as n
import glob

class LightconeRemap :
	__zMean = 0.88
	__Lbox = 2500.0 * uu.Mpc # Mpc/h
	__wdir = "/home2/jcomparat/eBOSS-LC/Multidark-lightcones/"
	__boxDir = "MD_2.5Gpc/"
	__snl =  ""
	__zsl =  None
	__zArray = n.arange(0.2,2.4,1e-1)
	__Hbox = 67.77 * uu.km / (uu.s * uu.Mpc)
	__outputName = __wdir + __boxDir + "lc_square.txt"
	__Melement = 23593750000.0

	def __init__(self,zMean,Lbox,wdir,boxDir,snl,zsl,zArray,Hbox,outputName,Melement ):
		self.__Lbox = Lbox # box length
		self.__Hbox = Hbox # Hubble constant at redshift 0 in the box
		self.__zMean = zMean # mean redshift desired
		self.__wdir = wdir # working directory
		self.__boxDir = boxDir # directory of the box where the snapshots a stored
		self.__snl = snl # snapshot list
		self.__zsl = zsl # corresponding redshift list
		self.__zArray = zArray # redshift for the dC - z conversion
		self.__outputName = outputName
		self.__Melement = Melement # mass of one particle in the box

	def set_Melement(self,Melement):
		self.__Melement = Melement

	def get_Melement(self):
		return self.__Melement 

	def set_zMean(self,zMean):
		self.__zMean = zMean

	def get_zMean(self):
		return self.__zMean 

	def set_outputName(self,outputName):
		self.__outputName = outputName

	def get_outputName(self):
		return self.__outputName 

	def set_zArray(self,zArray):
		self.__zArray = zArray

	def get_zArray(self):
		return self.__zArray 

	def set_Lbox(self,Lbox):
		self.__Lbox = Lbox

	def get_Lbox(self):
		return self.__Lbox 

	def set_Hbox(self,Hbox):
		self.__Hbox = Hbox

	def get_Hbox(self):
		return self.__Hbox 

	def set_wdir(self,wdir):
		self.__wdir = wdir

	def get_wdir(self):
		return self.__wdir 

	def set_boxDir(self,boxDir):
		self.__boxDir = boxDir

	def get_boxDir(self):
		return self.__boxDir 

	def set_snl(self,snl):
		self.__snl = snl

	def get_snl(self):
		return self.__snl 

	def set_zsl(self,zsl):
		self.__zsl = zsl

	def get_zsl(self):
		return self.__zsl 

	def defineLClimits(self,massFactor=50) :
		self.D_bar = c2.comoving_distance( self.get_zsl()[n.searchsorted( self.get_zsl(), self.get_zMean() )] ) * c2.h # Mpc / h
		print "D_bar=",self.D_bar
		#self.theta_max = n.arcsin(1. / ( 1. + 2. * self.D_bar/ self.get_Lbox() ) ) # * 180. / n.pi
		#print "theta_max=",self.theta_max
		#self.tg_theta_max = n.tan(self.theta_max)
		#print "tan theta_max=",self.tg_theta_max
		self.D_min = ( (self.D_bar - self.get_Lbox()/4.)**2. + self.get_Lbox()**2.  ) ** 0.5
		print "D_min=",self.D_min
		self.D_max = self.D_bar + self.get_Lbox()/4.
		print "D_max=",self.D_max
		self.dc_to_z = interp1d(c2.comoving_distance(self.get_zArray()) * c2.h ,self.get_zArray())
		print "dc to z interpolated for",self.get_zArray()
		self.Mmin=massFactor*self.get_Melement()

	def defineSnapshotLimits(self):
		# defines the limits of the LC
		dsl = c2.comoving_distance(self.get_zsl()) * c2.h # Mpc / h
		print "dsl=", dsl
		selection = (dsl > self.D_min) & (dsl < self.D_max)
		self.set_snl(self.get_snl()[selection])
		self.set_zsl(self.get_zsl()[selection])
		self.dTr=n.array((dsl[selection][1:]+dsl[selection][:-1])/2.)*uu.Mpc
		print "dTr=",self.dTr
		#trs=n.hstack(( self.D_min, dTr, self.D_max )) 
		#print "trs", trs
		self.D_transition = n.empty_like(n.arange(len(self.dTr)+2)) * uu.Mpc
		self.D_transition[0] = self.D_min
		self.D_transition[-1] = self.D_max
		self.D_transition[1:-1] = self.dTr
		print "D_transition=",self.D_transition 
	
	def selectline(self,line,dmin,dmax):
		# from a line outputs transformed coordinates ra, dec, R
		# based on standard Rockstar halo catalog structure
		x0=float(line[17]) * uu.Mpc #-self.get_Lbox()/2.
		y0=float(line[18]) * uu.Mpc #-self.get_Lbox()/2.
		z0=float(line[19]) * uu.Mpc #-self.get_Lbox()/2.+self.D_bar
		if z>self.get_Lbox()/2. :
			x = x0 - self.get_Lbox()
			y = y0 - self.get_Lbox()/2.
			z = z0 - self.get_Lbox()*3./4. + self.D_bar 
		else:
			x = x0 
			y = y0 - self.get_Lbox()/2.
			z = z0 - self.get_Lbox()/4. + self.D_bar 

		mvir=float(line[10])
		rr=x**2 + y**2 + z**2
		print line
		print x,y,z,rr**0.5, dmin, dmax, "-----",mvir,self.Mmin #, abs(x/z), abs(y/z), self.tg_theta_max
		if rr >= dmin**2. and rr < dmax**2. and mvir > self.Mmin: # and abs(x/z) < self.tg_theta_max and abs(y/z) < self.tg_theta_max :
			return 1.
		else:
			return 0.

	def transformLine(self,line):
		# if the point is chosen, the line is transformed into the LC.
		x0=float(line[17]) * uu.Mpc #-self.get_Lbox()/2.
		y0=float(line[18]) * uu.Mpc #-self.get_Lbox()/2.
		z0=float(line[19]) * uu.Mpc #-self.get_Lbox()/2.+self.D_bar
		if z>self.get_Lbox()/2. :
			x = x0 - self.get_Lbox()
			y = y0 - self.get_Lbox()/2.
			z = z0 - self.get_Lbox()*3./4. + self.D_bar 
		else:
			x = x0 
			y = y0 - self.get_Lbox()/2.
			z = z0 - self.get_Lbox()/4. + self.D_bar 

		a=float(line[0]) # scale factor
		ra=n.arctan(x/z)*180/n.pi
		dec=n.arctan(y/z)*180/n.pi

		vx=float(line[20]) * uu.km / uu.s
		vy=float(line[21]) * uu.km / uu.s
		vz=float(line[22]) * uu.km / uu.s

		rr=(x**2 + y**2 + z**2)**0.5
		redshiftR = self.dc_to_z(rr)

		vPara = (vx * x + vy * y + vz * z )/rr 
		rs = rr + vPara / (a * self.get_Hbox())
		redshiftS = self.dc_to_z(rs)
		#print a,x,y,z,ra,dec,vx,vy,vz,rr,vPara, self.get_Hbox(),redshiftR,redshiftS
		
		indexes = [1,4,5,6,10,11,12,13,16,26,37,38,44,45,57,58,59,60,61,63,64,67]
		return n.hstack((str(ra.value),str(dec.value),str(redshiftR),str(redshiftS),str(vPara.value),n.array(line)[indexes]))

	def constructLC(self):
		fo=open(self.get_outputName(),'a')
		for ii in range(len(self.get_snl())):
			fl=fileinput.input(self.get_snl()[ii])
			dmin=self.D_transition[ii]
			dmax=self.D_transition[ii+1]
			print "opens file ",self.get_snl()[ii],"and gets halo in ",dmin,dmax
			# self.tg_theta_max
			count=0
			for line in fl:
				if line[0]=="#" :
					count+=1
					continue

				yn=self.selectline(line.split(),dmin,dmax)
				if yn==0. :
					count+=1					
					continue

				if yn==1.:
					trL=self.transformLine(line.split())
					toWrite = " ".join(trL)+" \n"
					print "got line"
					print toWrite
					print count
					fo.write(toWrite)
					count+=1

			
			fl.close()

		fo.close()


wdir = "/home2/jcomparat/eBOSS-LC/Multidark-lightcones/"
boxDir = "MD_2.5Gpc/"
snlInter=n.array(glob.glob(wdir+boxDir+"h*.list"))
zInter=n.array([ 1/float(el.split('_')[-1][:-5])-1 for el in snlInter ])
ids=n.argsort(zInter)
snl=snlInter[ids]
zs=zInter[ids]
zArray=n.arange(0.1,2.4,1e-4)
Hbox=67.77 * uu.km / (uu.s * uu.Mpc)
outputName = wdir + boxDir + "MD_2.5Gpc_lc_remap_z0.88.txt"
Melement = 23593750000.0

lc1=LightconeRemap(zMean=0.88,Lbox=2500.0 * uu.Mpc,wdir=wdir,boxDir=boxDir,snl=snl,zsl=zs,zArray=zArray,Hbox=Hbox,outputName=outputName,Melement = Melement)
print "lc object created, defines limits now: "
lc1.defineLClimits(massFactor=50)
print "limits ok"
lc1.defineSnapshotLimits()
print "starts lc building"
lc1.constructLC()


# c2.luminosity_distance(1)*c2.h in Mpc / h
# c2.comoving_distance(1)*c2.h in Mpc / h

#c1=co.FlatLambdaCDM( H0=100.* uu.km / (uu.Mpc *uu.s), Om0=0.307, Tcmb0=2.725 *uu.K, Neff=3.05, m_nu=[ 0. ,   0. ,   0.06] *uu.eV, Ob0=0.0483)

#fileinput.input("1Gpc_3840_Planck1/BDM/CatshortV.0023.00.DAT")

#fl=fileinput.input("/data2/users/gustavo/BigMD/2.5_3840_Planck1/ROCKSTAR/out_29.list")

#name = "xxx"

#print(name)


