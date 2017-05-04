import glob
import os
import time
import numpy as n
import sys

def create_sat_files(fileName, outFileName ):
	#outFileName = fileName[:-5]+"_sat.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tpipe ifmt=fits in="""+fileName+""" cmd='select "pid>=0"' cmd='replacecol pid toInteger(pid)' cmd='replacecol id toInteger(id)' omode=out ofmt=fits out="""+outFileName
	os.system(command)

def create_cen_files(fileName, outFileName ):
	#outFileName = fileName[:-5]+"_cen.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tpipe ifmt=fits in="""+fileName+""" cmd='select "pid==-1 && mvir>13.5"' cmd='delcols "pid"' cmd='replacecol id toInteger(id)'  omode=out ofmt=fits out="""+outFileName
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
	satFileName = fileName[:-5] + "_sat_all.fits"
	cenFileName = fileName[:-5] + "_cen_all.fits"
	sat_in_sat_file = fileName[:-5] + "_subhalos_inSat.fits"
	outFileNameB = fileName[:-5] + "_subhalos_inSat2.fits"
	outFileNameA = fileName[:-5] + "_subhalos_inDistinct2.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+cenFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_1" suffix2="_cen" values1=pid_sat_n_1 values2=id omode=out ofmt=fits out="""+outFileNameA
	os.system(command)
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+satFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_1" suffix2="_sat_n_2" values1=pid_sat_n_1 values2=id omode=out ofmt=fits out="""+outFileNameB
	os.system(command)


def match_sat_cen_d3(fileName):
	# deprecated
	satFileName = fileName[:-5] + "_sat_all.fits"
	cenFileName = fileName[:-5] + "_cen_all.fits"
	sat_in_sat_file = fileName[:-5] + "_subhalos_inSat2.fits"
	outFileNameB = fileName[:-5] + "_subhalos_inSat3.fits"
	outFileNameA = fileName[:-5] + "_subhalos_inDistinct3.fits"
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+cenFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_2" suffix2="_cen" values1=pid_sat_n_2 values2=id omode=out ofmt=fits out="""+outFileNameA
	os.system(command)
	command = """java -jar /home/comparat/software/linux/stilts/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+satFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_2" suffix2="_sat_n_3" values1=pid_sat_n_2 values2=id omode=out ofmt=fits out="""+outFileNameB
	os.system(command)

def process_MD(files, outs):
	t0=time.time()
	for file in files:
		print file
		create_sat_files(file)
		print "create sat", time.time()-t0
		create_cen_files(file)
		print "create cen", time.time()-t0

	print "-----------------------------------------"	
	t0=time.time()
	for file in outs:
		print file
		concat_sat_files(file)
		print "concat sat", time.time()-t0
		concat_cen_files(file)
		print "concat cen", time.time()-t0
		
	print "-----------------------------------------"	
	t0=time.time()
	for file in outs:
		match_sat_cen(file)
		print "match", time.time()-t0
		

def match_cats(outs):
	t0=time.time()
	for file in outs:
		match_sat_cen(file)
		print "match", time.time()-t0
		

env = os.environ['MD10']

files = n.array(glob.glob(os.path.join(env,"snapshots", "out_0.74980_SAM_Nb_*.fits")))
files.sort()
"""
#first separate distinct and subhalos
for fil in files:
	bn = os.path.basename(fil)
	out_s = os.path.join(env,"substructure", bn[:-5]+"_subH.fits")
	out_d = os.path.join(env,"substructure", bn[:-5]+"_disT.fits")
	create_sat_files(fil, out_s)
	create_cen_files(fil, out_d)
"""
# concatenates the sat and central files into single files
fileList = os.path.join(env,"substructure", "*_disT.fits")
outFileName =  os.path.join(env,"substructure", "out_0.74980_disT_mvir_gt_135.fits")
concat_cen_files(fileList, outFileName)
"""
fileList = os.path.join(env,"substructure", "*_subH.fits")
outFileName =  os.path.join(env,"substructure", "out_0.74980_subH.fits")
concat_sat_files(fileList, outFileName)
"""
satFileName = os.path.join(env,"substructure", "out_0.74980_subH.fits")
cenFileName =  os.path.join(env,"substructure", "out_0.74980_disT_mvir.gt.135.fits")

match_sat_cen(satFileName, cenFileName, outFileNameA="out_0.74980_subH_inDistinct_d1.fits", outFileNameB="out_0.74980_subH_inSub_d1.fits")

#match_sat_cen(outs[-1])
#match_sat_cen_d2(outs[-1])
#match_sat_cen_d3(outs[-1])

