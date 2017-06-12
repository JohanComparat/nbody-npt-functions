from selectHaloCatalogsLib25 import *

# UNDER CONSTRUCTION
"""
ids=[]
for ii in range(len(zmin)):
	nG=int(nGal_Deg2[ii]* area)+1
	print "gets all halos for ",zmin[ii],"<z<",zmax[ii], "with mvir to mock ", nG, " galaxies."	
	IDhz,QTY,nn,bb=get_distrib_QTY(hdu, colN='mvir', zmin=zmin[ii], zmax=zmax[ii])
	ids.append(sham(nG,IDhz, QTY, nn,bb))

writerCatsAll("SHAM_complete_ELG_240.cat",n.sort(n.hstack((ids))))
"""

icPar=2.
ids=[]
for ii in range(len(zmin)):
        nG=int(nGal_Deg2[ii]* area)+1
        print "gets all halos for ",zmin[ii],"<z<",zmax[ii], "with mvir to mock ", nG, " galaxies."
        IDhz,QTY,nn,bb=get_distrib_QTY(hdu, colN='mvir', zmin=zmin[ii], zmax=zmax[ii])
	ids.append(shamIncomplete(icPar, nG,IDhz, QTY, nn,bb))

writerCatsAll("SHAM_incomplete_"+str(icPar)+"_ELG_240.cat",n.sort(n.hstack((ids))))


icPar=1.5
ids=[]
for ii in range(len(zmin)):
        nG=int(nGal_Deg2[ii]* area)+1
        print "gets all halos for ",zmin[ii],"<z<",zmax[ii], "with mvir to mock ", nG, " galaxies."
        IDhz,QTY,nn,bb=get_distrib_QTY(hdu, colN='mvir', zmin=zmin[ii], zmax=zmax[ii])
        ids.append(shamIncomplete(icPar, nG,IDhz, QTY, nn,bb))

writerCatsAll("SHAM_incomplete_"+str(icPar)+"_ELG_240.cat",n.sort(n.hstack((ids))))


import sys
sys.exit()



sats=n.arange(0.,0.31,0.025)
mMean=n.array([4,8,16,32,64,128,256])*10**(12)
ps=[]
ids,oN=[],[]
for el in mMean:
	for fs in sats:
		ps.append([el,el,fs])
		ids.append([])
		oN.append("")
		ps.append([el,el/2.,fs])
		ids.append([])
		oN.append("")
		ps.append([el,el/4.,fs])
		ids.append([])
		oN.append("")

print len(ps), len(ids), len(oN)

for ii in range(len(zmin)):
	nG=int(nGal_Deg2[ii]* area)+1
	print "gets all halos for ",zmin[ii],"<z<",zmax[ii], "with col5 to mock ", nG, " galaxies." 
	IDhz_c,QTY_c,nn_c,bb_c=get_distrib_QTY_cen(hdu, colN='mvir', zmin=zmin[ii], zmax=zmax[ii])
	IDhz_s,QTY_s,nn_s,bb_s=get_distrib_QTY_sat(hdu, colN='mvir', zmin=zmin[ii], zmax=zmax[ii])
	for jj in range(len(ps)):
		p1,p2,p3=ps[jj] #4e13,1e13,0.5111
		oN[jj]="SHAM_norm-mean"+str(p1)+"-sig"+str(p2)+"-fsat"+str(p3)+"_ELGdensity"
		ids[jj].append(selectGaussian_fsat(p1,p2,p3, nG,IDhz_c, QTY_c,IDhz_s, QTY_s ))

for jj in range(len(ps)):
	writerCatsAll(oN[jj],n.sort(n.hstack((ids[jj]))))

"""
idm1,idm2,idm3,idm4,idm5,idm6,idm7,idm8=[],[],[],[],[],[],[],[]
for ii in range(len(zmin)):
	nG=int(nGal_Deg2[ii]* area)+1
	print "gets all halos for ",zmin[ii],"<z<",zmax[ii], "with col5 to mock ", nG, " galaxies." 
	IDhz,QTY,nn,bb=get_distrib_QTY(hdu, colN='mvir', zmin=zmin[ii], zmax=zmax[ii])
	p1,p2=4e12,1e12
	idm1.append(selectGaussian(p1,p2, nG,IDhz, QTY, nn,bb))
	idm2.append(selectGaussian(p1,p2/2., nG,IDhz, QTY, nn,bb))

	p3,p4=2e12,5e11
	idm3.append(selectGaussian(p3,p4, nG,IDhz, QTY, nn,bb))
	idm4.append(selectGaussian(p3,p4/2., nG,IDhz, QTY, nn,bb))

	p5,p6=1e12,2.5e11
	idm5.append(selectGaussian(p5,p6, nG,IDhz, QTY, nn,bb))
	idm6.append(selectGaussian(p5,p6/2., nG,IDhz, QTY, nn,bb))

	p7,p8=8e11,2e11
	idm7.append(selectGaussian(p7,p8, nG,IDhz, QTY, nn,bb))
	idm8.append(selectGaussian(p7,p8/2., nG,IDhz, QTY, nn,bb))

writerCatsAll("SHAM_norm"+str(p1)+"sig"+str(p2)+"_ELGdensity",n.sort(n.hstack((idm1))))
writerCatsAll("SHAM_lognorm"+str(p1)+"sig"+str(p2)+"o2_ELGdensity",n.sort(n.hstack((idm2))))

writerCatsAll("SHAM_norm"+str(p3)+"sig"+str(p4)+"_ELGdensity",n.sort(n.hstack((idm3))))
writerCatsAll("SHAM_lognorm"+str(p3)+"sig"+str(p4)+"o2_ELGdensity",n.sort(n.hstack((idm4))))

writerCatsAll("SHAM_norm"+str(p5)+"sig"+str(p6)+"_ELGdensity",n.sort(n.hstack((idm5))))
writerCatsAll("SHAM_lognorm"+str(p5)+"sig"+str(p6)+"o2_ELGdensity",n.sort(n.hstack((idm6))))

writerCatsAll("SHAM_norm"+str(p7)+"sig"+str(p8)+"_ELGdensity",n.sort(n.hstack((idm7))))
writerCatsAll("SHAM_lognorm"+str(p7)+"sig"+str(p8)+"o2_ELGdensity",n.sort(n.hstack((idm8))))
"""






#print " does sham "
#id_mock_halos.append(sham( nGal=nG, IDhz=IDhz, QTY=QTY, nn=nn,bb=bb))
"""
print " does sham inc2"
id_mock_halos_inc2.append(shamIncomplete(incompFactor=10000., nGal=nG, IDhz=IDhz, QTY=QTY, nn=nn,bb=bb))
print " does sham inc5"
id_mock_halos_inc5.append(shamIncomplete(incompFactor=7000.,nGal=nG, IDhz=IDhz, QTY=QTY, nn=nn,bb=bb))
print " does sham inc10"
id_mock_halos_inc10.append(shamIncomplete(incompFactor=6000.,nGal=nG, IDhz=IDhz, QTY=QTY, nn=nn,bb=bb))
print " does sham inc100"
id_mock_halos_inc100.append(shamIncomplete(incompFactor=4000.,nGal=nG, IDhz=IDhz, QTY=QTY, nn=nn,bb=bb))
"""

"""
#writerCatsAll("SHAM_ELGdensity",n.sort(n.hstack((id_mock_halos))))
writerCatsAll("SHAM_inc10000_ELGdensity",n.sort(n.hstack((id_mock_halos_inc2))))
writerCatsAll("SHAM_inc7000_ELGdensity",n.sort(n.hstack((id_mock_halos_inc5))))
writerCatsAll("SHAM_inc6000_ELGdensity",n.sort(n.hstack((id_mock_halos_inc10))))
writerCatsAll("SHAM_inc4000_ELGdensity",n.sort(n.hstack((id_mock_halos_inc100))))

id_mock_halos.append(selectGaussian(1e12,5e11, nG,IDhz, QTY, nn,bb))

writerCatsAll("Halo_norm_selection_ELGdensity",n.sort(n.hstack((id_mock_halos))))

"""


