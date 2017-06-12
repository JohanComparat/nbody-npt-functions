#!/bin/bash/

#wdir = "/home2/jcomparat/eBOSS-LC/Multidark-lightcones/"
#boxDir = "MD_2.5Gpc/lc_square_0.6z1.2/"
#SNlist = glob.glob.(wdir + boxDir + "shell_*")

cd /home2/jcomparat/eBOSS-LC/Multidark-lightcones/MD_2.5Gpc/lc_square_0.6z1.2/

for i in shell_*.list; do
    java -jar ~/code/stilts.jar tpipe in=$i ifmt=ascii cmd='select "abs(ra)<30 && abs(dec)<30 && z>0.6 && z<1.2"' omode=out ofmt=fits out=$i.fits 
done

ls shell_*.fits > shell_list
java -jar ~/code/stilts.jar tcat in=@shell_list ifmt=fits omode=out ofmt=fits out=lightcone_MD_2.5Gpc_0.6z1.2.fits 
