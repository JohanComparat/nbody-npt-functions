import random
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import astropy.io.fits as pf
import numpy as n
from scipy.interpolate import interp1d
import scipy.stats as st

mockDir="/data2/DATA/eBOSS/ELG/HAM-mocks-decam240/"
NZfile="/home2/jcomparat/pyMDMC/NZ/decam240-radeczw.data.NZ"
# Cat_W134.cat.hNZ"
# "/home2/jcomparat/pyMDMC/NZ/vipers-W14-radeczW-06z10.cat.masked.NZ"
# "/home2/jcomparat/pyMDMC/NZ/Cat_W134.cat.hNZ"
print "outputs stored : ",mockDir
print "N(z) from ", NZfile

lcDir="/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc_new_rockS/"
lcFile=lcDir+"eboss-lc-newMD-z0.6-1.0.fits"
# loads the LC
hdu=pf.open(lcFile)
print lcFile, " loaded, columns:"
print hdu[1].data.dtype
area=(2*11.6929)**2.
#  ra(1) dec(2) zobs(3) zreal(4) Mvir(5) M200b(6) 
#   M200c(7) Macc(8) Mpeak(9) Vmax(10) Vacc(11) 
#   Vpeak(12) Vpec(13) Rs(14) ID(15) PID(16) 
#   UPID(17) SnapNUM(18)

# loads the randoms
# raR,decR=n.loadtxt(lcDir+"eboss-lc-newMD-z0.6-1.0.ran.cat",unpack=True)
cen=(hdu[1].data['col16']==-1)
sat=(cen==False)

sZ=lambda hdu,minz,maxz : (hdu[1].data['col3']>=minz)&(hdu[1].data['col3']<maxz)


zminIN,zmaxIN,nGal_Deg2IN,nGal_Mpc3IN=n.loadtxt(NZfile,unpack=True,usecols=(0,1,2,3))
ok=(nGal_Deg2IN>0)
zmin,zmax,nGal_Deg2,nGal_Mpc3 = zminIN[ok],zmaxIN[ok],nGal_Deg2IN[ok],nGal_Mpc3IN[ok]

def writerCats(name,idSel):
        print "writes ", name
        n.savetxt(mockDir+name+"_radecz.cat",n.transpose([hdu[1].data['col1'][idSel], hdu[1].data['col2'][idSel], hdu[1].data['col3'][idSel]]),fmt='%.8f %.8f %.5f')

def writerCatsAll(name,idSel):
        print "writes ", name
        n.savetxt(mockDir+name+"_radecz.cat",n.transpose([hdu[1].data['col1'][idSel], hdu[1].data['col2'][idSel], hdu[1].data['col3'][idSel]]))
        n.savetxt(mockDir+name+"_allCols.cat.gz" ,n.transpose([hdu[1].data['col1'][idSel], hdu[1].data['col2'][idSel], hdu[1].data['col3'][idSel], hdu[1].data['col5'][idSel], hdu[1].data['col16'][idSel]]))

def get_distrib_QTY(hdu, colN, zmin, zmax):
	IDh=n.arange(len(hdu[1].data[colN]))
	zsel=sZ(hdu,zmin,zmax)
	IDhz=IDh[zsel] # all ids in this redshift bin
	QTY=hdu[1].data[colN][zsel] # all QTY in this redshift bin
	nn,bb,pp=p.hist(QTY,cumulative=True,bins=len(QTY)/100)
	print zmin,zmax,len(IDhz)
	return IDhz,QTY,nn,bb

def get_distrib_QTY_cen(hdu, colN, zmin, zmax):
	IDh=n.arange(len(hdu[1].data[colN]))
	zsel=sZ(hdu,zmin,zmax)&(cen)
	IDhz=IDh[zsel] # all ids in this redshift bin
	QTY=hdu[1].data[colN][zsel] # all QTY in this redshift bin
	nn,bb,pp=p.hist(QTY,cumulative=True,bins=len(QTY)/100)
	print zmin,zmax,len(IDhz)
	return IDhz,QTY,nn,bb

def get_distrib_QTY_sat(hdu, colN, zmin, zmax):
	IDh=n.arange(len(hdu[1].data[colN]))
	zsel=sZ(hdu,zmin,zmax)&(sat)
	IDhz=IDh[zsel] # all ids in this redshift bin
	QTY=hdu[1].data[colN][zsel] # all QTY in this redshift bin
	nn,bb,pp=p.hist(QTY,cumulative=True,bins=len(QTY)/100)
	print zmin,zmax,len(IDhz)
	return IDhz,QTY,nn,bb

def sham(nGal,IDhz, QTY, nn,bb):
	mfc=interp1d(nn,(bb[1:]+bb[:-1])/2.)
	QTYmax=mfc(len(QTY))
	QTYmin=mfc(len(QTY)-nGal)
	qsel=(QTY>QTYmin)&(QTY<=QTYmax)
	IDhzq=IDhz[qsel]
	print zmin,zmax,nGal,len(IDhzq)
	return IDhzq

def shamIncomplete(incompFactor, nGal,IDhz, QTY, nn,bb):
	mfc=interp1d(nn,(bb[1:]+bb[:-1])/2.)
	mfcInv=interp1d((bb[1:]+bb[:-1])/2.,nn)
	QTYmaxAll=mfc(len(QTY))/incompFactor
	Nmax=mfcInv(QTYmaxAll)
	QTYmax=mfc(Nmax)
	QTYmin=mfc(Nmax-nGal)
	qsel=(QTY>QTYmin)&(QTY<=QTYmax)
	IDhzq=IDhz[qsel]
	return IDhzq

def sham_QTY_max(QTY_max, nGal,IDhz, QTY, nn,bb):
	mfc=interp1d(nn,(bb[1:]+bb[:-1])/2.)
	mfcInv=interp1d((bb[1:]+bb[:-1])/2.,nn)
	Nmax=mfcInv(QTY_max)
	QTYmax=mfc(Nmax)
	QTYmin=mfc(Nmax-nGal)
	qsel=(QTY>QTYmin)&(QTY<=QTYmax)
	IDhzq=IDhz[qsel]
	return IDhzq

def selectGaussian(position,scatter, nGal,IDhz, QTY, nn,bb):
	# constructs the QTY intervals around the distribution
	expected_cdf=lambda x : st.norm.cdf(x, loc=position, scale=scatter)
	interval = [ position - 9 * scatter , position + 9 * scatter]
	xs=n.arange(interval[0],interval[1],(interval[1]-interval[0])/1000.)
	out=expected_cdf(xs)
	expected_cdf_inv=interp1d(out,xs)
	boundaries=n.hstack((expected_cdf_inv(0.01),expected_cdf_inv(n.arange(0.1,0.91,0.1)), interval[1]))
	# gets the number of halos to select
	expected_cdf_tot=lambda x : nGal * st.norm.cdf(x, loc=position, scale=scatter)
	Up=expected_cdf_tot(boundaries[1:])
	Low=n.hstack(( 0., expected_cdf_tot(boundaries[1:])[:-1] ))
	N2select=Up-Low
	print N2select,Up,Low
	# select in mass in the box
	qsels=n.array([ (QTY>boundaries[ii])&(QTY<=boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
	IDhzqAll=n.array([ IDhz[qs] for qs in qsels ])
	# random downsample to the N2select in each bin
	i=0
	ids_selected=[]
	for arr in IDhzqAll:
		random.shuffle(arr)
		ids_selected.append(arr[:N2select[i]])
		i+=1

	ids_selected=n.hstack(( n.array(ids_selected) ))
	return ids_selected

def selectGaussian_fsat(position,scatter,fsat, nGal,IDhz_c, QTY_c,IDhz_s, QTY_s ):
	nSat=int(nGal*fsat)
	print nGal,nSat,fsat,position,scatter
	# constructs the QTY intervals around the distribution 
	expected_cdf=lambda x : st.norm.cdf(x, loc=position, scale=scatter)
	interval = [ position - 9 * scatter , position + 9 * scatter]
	xs=n.arange(interval[0],interval[1],(interval[1]-interval[0])/1000.)
	out=expected_cdf(xs)
	expected_cdf_inv=interp1d(out,xs)
	boundaries=n.hstack((expected_cdf_inv(0.01),expected_cdf_inv(n.arange(0.1,0.91,0.1)), interval[1]))
	# gets the number of halos to select the SAT
	print "sats"
	expected_cdf_s = lambda x : nSat * st.norm.cdf(x, loc=position, scale=scatter)
	Up_s = expected_cdf_s(boundaries[1:])
	Low_s = n.hstack(( 0., expected_cdf_s(boundaries[1:])[:-1] ))
	N2select_s = Up_s-Low_s
	# select in mass in the box
	qsels_s=n.array([ (QTY_s>boundaries[ii])&(QTY_s<=boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
	IDhzqAll_s=n.array([ IDhz_s[qs] for qs in qsels_s ])

	# random downsample to the N2select in each bin
	i=0
	ids_selected_s=[]
	for arr2 in IDhzqAll_s:
		random.shuffle(arr2)
		print len(arr2),int(N2select_s[i])
		ids_selected_s.append(arr2[:int(N2select_s[i])])
		i+=1

	id_s=n.hstack((n.array(ids_selected_s)))
	
	nSatReal=len(id_s)
	nCen=nGal-nSatReal
	print nGal,nSat,nCen,fsat,position,scatter

	# gets the number of halos to select the CEN, compatible with the sat fraction to get the right density.
	print "centrals"
	expected_cdf_c = lambda x : nCen * st.norm.cdf(x, loc=position, scale=scatter)
	Up_c = expected_cdf_c(boundaries[1:])
	Low_c = n.hstack(( 0., expected_cdf_c(boundaries[1:])[:-1] ))
	N2select_c = Up_c-Low_c
	# select in mass in the box
	qsels_c=n.array([ (QTY_c>boundaries[ii])&(QTY_c<=boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
	IDhzqAll_c=n.array([ IDhz_c[qs] for qs in qsels_c ])

	# random downsample to the N2select in each bin
	i=0
	ids_selected_c=[]
	for arr in IDhzqAll_c:
		random.shuffle(arr)
		print len(arr),int(N2select_c[i])
		ids_selected_c.append(arr[:int(N2select_c[i])])
		i+=1

	id_c=n.hstack((n.array(ids_selected_c)))

	print len(id_c),len(id_s)
	ids_selected=n.hstack((id_c,id_s ))
	print len(ids_selected)
	return ids_selected




def selectLogNorm(position,scatter, nGal,IDhz, QTY, nn,bb):
	# constructs the QTY intervals around the distribution
	expected_cdf=lambda x : st.lognorm.cdf(x, position, scatter)
	interval = [ position - 9 * scatter , position + 9 * scatter]
	xs=n.arange(interval[0],interval[1],(interval[1]-interval[0])/1000.)
	out=expected_cdf(xs)
	expected_cdf_inv=interp1d(out,xs)
	boundaries=n.hstack((expected_cdf_inv(0.01),expected_cdf_inv(n.arange(0.1,0.91,0.1)), interval[1]))
	# gets the number of halos to select
	expected_cdf_tot=lambda x : nGal * st.lognorm.cdf(x, position, scatter)
	Up=expected_cdf_tot(boundaries[1:])
	Low=n.hstack(( 0., expected_cdf_tot(boundaries[1:])[:-1] ))
	N2select=Up-Low
	print N2select,Up,Low
	# select in mass in the box
	qsels=n.array([ (QTY>boundaries[ii])&(QTY<=boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
	IDhzqAll=n.array([ IDhz[qs] for qs in qsels ])
	# random downsample to the N2select in each bin
	i=0
	ids_selected=[]
	for arr in IDhzqAll:
		random.shuffle(arr)
		ids_selected.append(arr[:N2select[i]])
		i+=1

	ids_selected=n.hstack(( n.array(ids_selected) ))
	return ids_selected

