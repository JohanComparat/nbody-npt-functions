import numpy as n
import os

def writeParamFile(dataPath,randomPath,type,decade=""):
    f=open(dataPath[:-4]+".param2PCF_"+type+decade,'a')
    f.write("data_filename= "+dataPath+" \n")
    f.write("random_filename= "+randomPath+" \n")
    f.write("input_format= 2 \n")
    f.write("mask_filename= 'none' \n")
    f.write("z_dist_filename= 'none' \n")
    f.write("output_filename= "+dataPath[:-4]+"_"+type+decade+".dat \n")
    f.write("num_lines= all \n")
    f.write("corr_type= "+type+" \n")
    f.write("corr_estimator= LS \n")
    f.write("np_rand_fact= 10 \n")
    f.write("omega_M= 0.307 \n")
    f.write("omega_L= 0.693 \n")
    f.write("w= -1 \n")
    f.write("radial_aperture= 1 \n")
    f.write("use_pm= 0 \n")
    f.write("n_pix_sph= 2048 \n")
    f.close()


deltaZ=0.1
zmin=0.6 
zmax=1.2
zs = n.arange(zmin, zmax+deltaZ/2., deltaZ)
finalDensity = n.ones_like(zs)[:-1]

def tomography(dataFile= "/data1/DATA/eBOSS/COMBINED/eBOSS-combined_LRG_QSO_ELG.dat", randomFile= "/data1/DATA/eBOSS/COMBINED/eBOSS-combined_LRG_QSO_ELG.ran", zs=zs , finalDensity=finalDensity ):
	"""
	cuts a file into multiple redshift bins.
	"""
	print "opens ", dataFile
	raD, decD, zD, typeD = n.loadtxt(dataFile, usecols=(0,1,3,10), unpack = True)
	print "opens ", randomFile
	raR, decR, zR, typeR = n.loadtxt(randomFile, usecols=(0,1,3,4) , unpack = True)
	print "selects "
	elg = (typeD == 3)
	elgR = (typeR == 3)
	rd_D = n.random.random(len(raD))
	for i in range(len(zs)-1):
		# now for the data
		sel =  (elg) & (zD>zs[i]) & (zD<zs[i+1]) & (rd_D <=finalDensity[i])
		data_path = dataFile[:-4]+"_"+str(n.round(zs[i],3))+"_z_"+ str(n.round(zs[i+1],3)) +".dat"
		print zs[i], zs[i+1], data_path
		n.savetxt(data_path,n.transpose([ raD[sel], decD[sel], zD[sel] ]))
		# now for the randoms
		sel =  (elgR) & (zR>zs[i]) & (zR<zs[i+1])
		random_path = dataFile[:-4]+"_"+str(n.round(zs[i],3))+"_z_"+ str(n.round(zs[i+1],3)) +".ran"
		print random_path
		n.savetxt(random_path,n.transpose([ raR[sel], decR[sel], zR[sel] ]))
		# writes clustering paramFiles
		writeParamFile(data_path, random_path, "monopole")
		writeParamFile(data_path, random_path, "angular", "_d1")
		writeParamFile(data_path, random_path, "angular", "_d2")
		writeParamFile(data_path, random_path, "angular", "_d3")


tomography()
