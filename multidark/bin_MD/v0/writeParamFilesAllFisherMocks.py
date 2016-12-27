import glob
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


mockDir="/data2/DATA/eBOSS/ELG/mocks_fischerGRIW1/"

haloList=glob.glob(mockDir+"mocks_fischerGRIW1SHAM*_radecz.cat")

os.system("rm -rf "+mockDir+"*param2PCF*")

for el in haloList:
	writeParamFile(el,el+".random.cat","monopole")
        writeParamFile(el,el+".random.cat","angular","_d1")
	writeParamFile(el,el+".random.cat","angular","_d2")
        writeParamFile(el,el+".random.cat","angular","_d3")

