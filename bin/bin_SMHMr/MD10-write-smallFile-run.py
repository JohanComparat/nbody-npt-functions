import sys
import glob
import os

for ii in range(46):
	print( "os.system('screen -S run_md_"+str(ii)+"')")

for ii in range(46):
	print( "python MD10-write-smallFile.py "+str(ii))

screen -S run_md_0
python MD10-write-smallFile.py 0
python MD10-write-smallFile.py 1
python MD10-write-smallFile.py 2
python MD10-write-smallFile.py 3
python MD10-write-smallFile.py 4
screen -S run_md_1
python MD10-write-smallFile.py 5
python MD10-write-smallFile.py 6
python MD10-write-smallFile.py 7
python MD10-write-smallFile.py 8
python MD10-write-smallFile.py 9
screen -S run_md_2
python MD10-write-smallFile.py 10
python MD10-write-smallFile.py 11
python MD10-write-smallFile.py 12
python MD10-write-smallFile.py 13
python MD10-write-smallFile.py 14
screen -S run_md_3
python MD10-write-smallFile.py 15
python MD10-write-smallFile.py 16
python MD10-write-smallFile.py 17
python MD10-write-smallFile.py 18
python MD10-write-smallFile.py 19
screen -S run_md_4
python MD10-write-smallFile.py 20
python MD10-write-smallFile.py 21
python MD10-write-smallFile.py 22
python MD10-write-smallFile.py 23
python MD10-write-smallFile.py 24
screen -S run_md_5
python MD10-write-smallFile.py 25
python MD10-write-smallFile.py 26
python MD10-write-smallFile.py 27
python MD10-write-smallFile.py 28
python MD10-write-smallFile.py 29
screen -S run_md_6
python MD10-write-smallFile.py 30
python MD10-write-smallFile.py 31
python MD10-write-smallFile.py 32
python MD10-write-smallFile.py 33
python MD10-write-smallFile.py 34
screen -S run_md_7
python MD10-write-smallFile.py 35
python MD10-write-smallFile.py 36
python MD10-write-smallFile.py 37
python MD10-write-smallFile.py 38
python MD10-write-smallFile.py 39
screen -S run_md_8
python MD10-write-smallFile.py 40
python MD10-write-smallFile.py 41
python MD10-write-smallFile.py 42
python MD10-write-smallFile.py 43
python MD10-write-smallFile.py 44
screen -S run_md_9
python MD10-write-smallFile.py 45


