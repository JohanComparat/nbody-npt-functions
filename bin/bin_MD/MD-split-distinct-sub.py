import glob
import os
import time
import numpy as n
import sys

def create_sat_files(fileName):
	outFileName = fileName[:-5]+"_sat.fits"
	command = """java -jar /home2/jcomparat/code/stilts.jar tpipe ifmt=fits in="""+fileName+""" cmd='select "pid>=0"' cmd='replacecol pid toInteger(pid)' cmd='replacecol id toInteger(id)' omode=out ofmt=fits out="""+outFileName
	os.system(command)

def create_cen_files(fileName):
	outFileName = fileName[:-5]+"_cen.fits"
	command = """java -jar /home2/jcomparat/code/stilts.jar tpipe ifmt=fits in="""+fileName+""" cmd='select "pid==-1"' cmd='delcols "pid x y z"' cmd='replacecol id toInteger(id)'  omode=out ofmt=fits out="""+outFileName
	os.system(command)

def concat_sat_files(fileName):
	os.system("ls "+fileName[:-5]+"*_sat.fits > list2Concat")
	outFileName = fileName[:-5] + "_sat_all.fits"
	command = """java -jar /home2/jcomparat/code/stilts.jar tcat ifmt=fits in=@list2Concat omode=out ofmt=fits out="""+outFileName
	os.system(command)
	os.system("rm list2Concat")
	
def concat_cen_files(fileName):
	os.system("ls "+fileName[:-5]+"*_cen.fits > list2Concat")
	outFileName = fileName[:-5] + "_cen_all.fits"
	command = """java -jar /home2/jcomparat/code/stilts.jar tcat ifmt=fits in=@list2Concat omode=out ofmt=fits out="""+outFileName
	os.system(command)
	os.system("rm list2Concat")
	
def match_sat_cen(fileName):
	satFileName = fileName[:-5] + "_sat_all.fits"
	cenFileName = fileName[:-5] + "_cen_all.fits"
	outFileNameA = fileName[:-5] + "_subhalos_inDistinct.fits"
	outFileNameB = fileName[:-5] + "_subhalos_inSat.fits"
	command = """java -jar /home2/jcomparat/code/stilts.jar tmatch2 ifmt1=fits in1="""+satFileName+""" ifmt2=fits in2="""+cenFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat" suffix2="_cen" values1=pid values2=id omode=out ofmt=fits out="""+outFileNameA
	os.system(command)
	command = """java -jar /home2/jcomparat/code/stilts.jar tmatch2 ifmt1=fits in1="""+satFileName+""" ifmt2=fits in2="""+satFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n" suffix2="_sat_n_1" values1=pid values2=id omode=out ofmt=fits out="""+outFileNameB
	os.system(command)

def match_sat_cen_d2(fileName):
	satFileName = fileName[:-5] + "_sat_all.fits"
	cenFileName = fileName[:-5] + "_cen_all.fits"
	sat_in_sat_file = fileName[:-5] + "_subhalos_inSat.fits"
	outFileNameB = fileName[:-5] + "_subhalos_inSat2.fits"
	outFileNameA = fileName[:-5] + "_subhalos_inDistinct2.fits"
	command = """java -jar /home2/jcomparat/code/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+cenFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_1" suffix2="_cen" values1=pid_sat_n_1 values2=id omode=out ofmt=fits out="""+outFileNameA
	os.system(command)
	command = """java -jar /home2/jcomparat/code/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+satFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_1" suffix2="_sat_n_2" values1=pid_sat_n_1 values2=id omode=out ofmt=fits out="""+outFileNameB
	os.system(command)


def match_sat_cen_d3(fileName):
	satFileName = fileName[:-5] + "_sat_all.fits"
	cenFileName = fileName[:-5] + "_cen_all.fits"
	sat_in_sat_file = fileName[:-5] + "_subhalos_inSat2.fits"
	outFileNameB = fileName[:-5] + "_subhalos_inSat3.fits"
	outFileNameA = fileName[:-5] + "_subhalos_inDistinct3.fits"
	command = """java -jar /home2/jcomparat/code/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+cenFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_2" suffix2="_cen" values1=pid_sat_n_2 values2=id omode=out ofmt=fits out="""+outFileNameA
	os.system(command)
	command = """java -jar /home2/jcomparat/code/stilts.jar tmatch2 ifmt1=fits in1="""+sat_in_sat_file+""" ifmt2=fits in2="""+satFileName+""" find=all matcher=exact join=1and2 fixcols=all suffix1="_sat_n_2" suffix2="_sat_n_3" values1=pid_sat_n_2 values2=id omode=out ofmt=fits out="""+outFileNameB
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
		match_sat_cen_d2(file)
		match_sat_cen_d3(file)
		print "match", time.time()-t0
		
files = n.array(glob.glob("/data2/DATA/eBOSS/DarkSkies/snapshots/ds14_*_PM_Nb_*.fits"))
files.sort()
outs = n.array(glob.glob("/data2/DATA/eBOSS/DarkSkies/snapshots/ds14_*.dat"))
outs.sort()
process_MD(files, outs)
match_cats(outs)
sys.exit()

files = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_0.4Gpc/snapshots/out_*_PM_Nb_?.fits"))
files.sort()
outs = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_0.4Gpc/snapshots/out_*.list"))
outs.sort()
#process_MD(files, outs)
match_cats(outs)

files = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc/snapshots/out_*_PM_Nb_?.fits"))
files.sort()
outs = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_1Gpc/snapshots/out_*.list"))
outs.sort()
#process_MD(files, outs)
match_cats(outs)


files = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/out_*_PM_Nb_?.fits"))
files.sort()
outs = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/out_*.list"))
outs.sort()
#process_MD(files, outs)
match_cats(outs)

files = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5GpcNW/snapshots/out_*_PM_Nb_?.fits"))
files.sort()
outs = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5GpcNW/snapshots/out_*.list"))
outs.sort()
#process_MD(files, outs)
match_cats(outs)

files = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_4GpcNW/snapshots/out_*_PM_Nb_?.fits"))
files.sort()
outs = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_4GpcNW/snapshots/out_*.list"))
outs.sort()
#process_MD(files, outs)
match_cats(outs)

files = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/out_*_PM_Nb_?.fits"))
files.sort()
outs = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_4Gpc/snapshots/out_*.list"))
outs.sort()
#process_MD(files, outs)
match_cats(outs)

