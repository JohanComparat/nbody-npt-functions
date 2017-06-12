from selectHaloCatalogsLib import *


ps=[10**14,10**13.9,10**13.8, 10**13.7,10**13.6,10**13.5, 10**13.4,10**13.3,10**13.2]
ids=[[],[],[], [],[],[], [],[],[]]
oN=["","","", "","","", "","",""]
for ii in range(len(zmin)):
	nG=int(nGal_Deg2[ii]* area)+1
	print "gets all halos for ",zmin[ii],"<z<",zmax[ii], "with col5 to mock ", nG, " galaxies." 
	IDhz,QTY,nn,bb=get_distrib_QTY(hdu, colN='col5', zmin=zmin[ii], zmax=zmax[ii])
	for jj in range(len(ps)):
		p1=ps[jj] #4e13,1e13,0.5111
		oN[jj]="SHAM_massMax-"+str(n.round(n.log10(p1),4))+"_VIPERSnz"
		ids[jj].append(sham_QTY_max(p1,nGal=nG, IDhz=IDhz, QTY=QTY, nn=nn,bb=bb))


for jj in range(len(ps)):
	writerCatsAll(oN[jj],n.sort(n.hstack((ids[jj]))))

