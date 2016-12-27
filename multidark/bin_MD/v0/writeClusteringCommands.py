import glob
import os
mockDir="/data2/DATA/eBOSS/ELG/HAM-mocks-decam240/"

# for the angular d2
os.system("rm -rf runAngd1.sh")
pl=glob.glob(mockDir+"*param*angular*d1")

f=open("runAngd1.sh",'a')
f.write("""#!/bin/bash/ \n""")
for el in pl:
    f.write("""/home2/jcomparat/code/CUTE-1.1A1/CUTE/CUTE """+el+"\n")


f.close()

# for the angular d2
os.system("rm -rf runAngd2.sh")
pl=glob.glob(mockDir+"*param*angular*d2")

f=open("runAngd2.sh",'a')
f.write("""#!/bin/bash/ \n""")
for el in pl:
    f.write("""/home2/jcomparat/code/CUTE-1.1A2/CUTE/CUTE """+el+"\n")


f.close()


# for the angular d3
os.system("rm -rf runAngd3.sh")
pl=glob.glob(mockDir+"*param*angular*d3")

f=open("runAngd3.sh",'a')
f.write("""#!/bin/bash/ \n""")
for el in pl:
    f.write("""/home2/jcomparat/code/CUTE-1.1A3/CUTE/CUTE """+el+"\n")


f.close()

# for the monopole
os.system("rm -rf runMono.sh")
pl=glob.glob(mockDir+"*param*monopole")

f=open("runMono.sh",'a')
f.write("""#!/bin/bash/ \n""")
for el in pl:
    f.write("""/home2/jcomparat/code/CUTE-1.1M/CUTE/CUTE """+el+"\n")


f.close()

