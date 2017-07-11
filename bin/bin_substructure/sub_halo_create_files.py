import glob
import os
import time
import numpy as n
import sys

def create_sat_files(fileName, outFileName ):
	#outFileName = fileName[:-5]+"_sat.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tpipe ifmt=fits in="""+fileName+""" cmd='select "x>0 && pid>=0"' omode=out ofmt=fits out="""+outFileName
	os.system(command)

def create_cen_files(fileName, outFileName ):
	#outFileName = fileName[:-5]+"_cen.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tpipe ifmt=fits in="""+fileName+""" cmd='select "x>0 && pid==-1 && mvir>13.5"' omode=out ofmt=fits out="""+outFileName
	os.system(command)

def concat_sat_files(fileList, outFileName):
	os.system("ls "+fileList+" > list2Concat")
	#outFileName = fileName[:-5] + "_sat_all.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tcat ifmt=fits in=@list2Concat omode=out ofmt=fits out="""+outFileName
	os.system(command)
	os.system("rm list2Concat")
	
def concat_cen_files(fileList, outFileName):
	os.system("ls "+fileList+" > list2Concat")
	#outFileName = fileName[:-5] + "_cen_all.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tcat ifmt=fits in=@list2Concat omode=out ofmt=fits out="""+outFileName
	os.system(command)
	os.system("rm list2Concat")
	
def match_sat_cen(satFileName, cenFileName, outFileNameA="_inDistinct_d1.fits", outFileNameB="_inSub_d1.fits"):
	#satFileName = fileName[:-5] + "_sat_all.fits"
	#cenFileName = fileName[:-5] + "_cen_all.fits"
	#outFileNameA = fileName[:-5] + "_subhalos_inDistinct.fits"
	#outFileNameB = fileName[:-5] + "_subhalos_inSat.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+satFileName+""" ifmt2=fits in2="""+cenFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat" suffix2="_cen" values1=pid values2=id omode=out ofmt=fits out="""+outFileNameA
	os.system(command)
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+satFileName+""" ifmt2=fits in2="""+satFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n" suffix2="_sat_n_1" values1=pid values2=id omode=out ofmt=fits out="""+outFileNameB
	os.system(command)

def match_sat_cen_d2(satFileName, cenFileName, sat_in_sat_file, outFileNameA="inDistinct_d2.fits", outFileNameB="inSub_d2.fits"):
	#satFileName = fileName[:-5] + "_sat_all.fits"
	#cenFileName = fileName[:-5] + "_cen_all.fits"
	#sat_in_sat_file = fileName[:-5] + "_subhalos_inSat.fits"
	#outFileNameB = fileName[:-5] + "_subhalos_inSat2.fits"
	#outFileNameA = fileName[:-5] + "_subhalos_inDistinct2.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+cenFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_1" suffix2="_cen" values1=pid_sat_n_1 values2=id omode=out ofmt=fits out="""+outFileNameA
	os.system(command)
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+satFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_1" suffix2="_sat_n_2" values1=pid_sat_n_1 values2=id omode=out ofmt=fits out="""+outFileNameB
	os.system(command)


def match_sat_cen_d3(satFileName, cenFileName, sat_in_sat_file, outFileNameA="inDistinct_d2.fits", outFileNameB="inSub_d2.fits"):
	# deprecated
	#satFileName = fileName[:-5] + "_sat_all.fits"
	#cenFileName = fileName[:-5] + "_cen_all.fits"
	#sat_in_sat_file = fileName[:-5] + "_subhalos_inSat2.fits"
	#outFileNameB = fileName[:-5] + "_subhalos_inSat3.fits"
	#outFileNameA = fileName[:-5] + "_subhalos_inDistinct3.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+cenFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_2" suffix2="_cen" values1=pid_sat_n_2 values2=id omode=out ofmt=fits out="""+outFileNameA
	os.system(command)
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+satFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_2" suffix2="_sat_n_3" values1=pid_sat_n_2 values2=id omode=out ofmt=fits out="""+outFileNameB
	os.system(command)


def process_MD(snap_num = "128p"):
	t0=time.time()
	env = os.environ['MD10']
	files = n.array(glob.glob(os.path.join(env,"work_agn", "out_"+snap_num+"_SAM_Nb_?.fits")))
	files.sort()
	if len(files)==0:
		print "no input files"
		return 0.
	else :
		print "start"
		#first separate distinct and subhalos
		for fil in files:
			print fil
			bn = os.path.basename(fil)
			out_s = os.path.join(env,"substructure", bn[:-5]+"_subH.fits")
			out_d = os.path.join(env,"substructure", bn[:-5]+"_disT.fits")
			create_sat_files(fil, out_s)
			print "sat", time.time()-t0
			create_cen_files(fil, out_d)
			print "cen", time.time()-t0

		
		# concatenates the sat and central files into single files
		fileList = os.path.join(env,"substructure", "*_disT.fits")
		cenFileName =  os.path.join(env,"substructure", "out_"+snap_num+"_disT_mvir_gt_135.fits")
		concat_cen_files(fileList, cenFileName)
		print "cat", time.time()-t0

		fileList = os.path.join(env,"substructure", "*_subH.fits")
		satFileName =  os.path.join(env,"substructure", "out_"+snap_num+"_subH_all.fits")
		concat_sat_files(fileList, satFileName)
		print "cat", time.time()-t0

		# 1st level match
		sat_in_cen_d1 =os.path.join(env, "substructure", "out_"+snap_num+"_subH_inDistinct_d1.fits")
		sat_in_sat_d1 =os.path.join(env, "substructure", "out_"+snap_num+"_subH_inSub_d1.fits")
		match_sat_cen(satFileName, cenFileName, sat_in_cen_d1 , sat_in_sat_d1)
		print "match", time.time()-t0
		
		# 2nd level match
		sat_in_cen_d2 = os.path.join(env, "substructure", "out_"+snap_num+"_subH_inDistinct_d2.fits")
		sat_in_sat_d2 = os.path.join(env, "substructure", "out_"+snap_num+"_subH_inSub_d2.fits")
		match_sat_cen_d2(satFileName, cenFileName, sat_in_sat_d1, sat_in_cen_d2, sat_in_sat_d2)
		print "match", time.time()-t0
		
		# 3rd level match
		sat_in_cen_d3 = os.path.join(env, "substructure", "out_"+snap_num+"_subH_inDistinct_d3.fits")
		sat_in_sat_d3 = os.path.join(env, "substructure", "out_"+snap_num+"_subH_inSub_d3.fits")
		match_sat_cen_d3(satFileName, cenFileName, sat_in_sat_d2, sat_in_cen_d3, sat_in_sat_d3)
		print "match", time.time()-t0
		
		return 1.
		
import astropy.io.fits as fits

sf = fits.open(os.path.join(os.environ["MD10"], "output_MD_1.0Gpc.fits"))[1].data

for snap in sf['snap_name']:
	print snap
	if os.path.isfile(os.path.join(os.environ["MD10"], "substructure", "out_"+snap+"_subH_inSub_d3.fits")) == False : 
		done = process_MD(snap)
	

os.system("rm out_*_SAM_Nb_*_disT.fits")
os.system("rm out_*_SAM_Nb_*_subH.fits")


#import astropy.io.fits as fits
#import numpy as n

#da = fits.open(sat_in_cen_d1)[1].data

#dist = ((da['x_sat']-da['x_cen'])**2. + (da['y_sat']-da['y_cen'])**2. + (da['z_sat']-da['z_cen'])**2.)**0.5

#inside =(dist < da['rvir_cen']/1000.)
