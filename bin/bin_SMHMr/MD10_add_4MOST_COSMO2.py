from populateHaloCatalogsLib import *

ps=[[1e12,5e11,0.1],[1e12,5e11,0.2],[1e12,5e11,0.3],[1e12,5e11,0.4],[1e12,5e11,0.5],[1e12,5e11,0.6]]
ids=[[],[],[], [],[],[]]
oN=["","","", "","",""]
for ii in range(len(zmin)):
	nG=int(nGal_Deg2[ii]* area)+1
	print "gets all halos for ",zmin[ii],"<z<",zmax[ii], "with col5 to mock ", nG, " galaxies." 
	IDhz_c,QTY_c,nn_c,bb_c=get_distrib_QTY_cen(hdu, colN='col5', zmin=zmin[ii], zmax=zmax[ii])
	IDhz_s,QTY_s,nn_s,bb_s=get_distrib_QTY_sat(hdu, colN='col5', zmin=zmin[ii], zmax=zmax[ii])
	for jj in range(len(ps)):
		p1,p2,p3=ps[jj] #4e13,1e13,0.5111
		oN[jj]="SHAM_norm-mean"+str(p1)+"-sig"+str(p2)+"-fsat"+str(p3)+"_ELGdensity"
		ids[jj].append(selectGaussian_fsat(p1,p2,p3, nG,IDhz_c, QTY_c,IDhz_s, QTY_s ))

for jj in range(len(ps)):
	writerCatsAll(oN[jj],n.sort(n.hstack((ids[jj]))))
