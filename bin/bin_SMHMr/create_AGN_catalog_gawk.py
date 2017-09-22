# overall python packages
import glob
#import astropy.io.fits as fits
# 2397897 143.540054 0.032711 20.449619 119.370173 9.753314 33.197590 -1.000000 25.191960 40.977921 2 127
# ------  --------   --------  ra        dec 

import os
import time
import numpy as n
import sys
t0=time.time()
#from astropy.cosmology import FlatLambdaCDM
#import astropy.units as u
#cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)


def get_AGN_catalog(env='MD10'):
        # gets the file list to add the Xray luminosity
        fileList = n.array(glob.glob(os.path.join(os.environ[env], "light-cone", "MDPL2_ROCKSTAR_FluxProj_*_000_AGN.dat" )))
        fileList.sort()
        #print fileList
        #print fileList[0]
        #data = n.loadtxt(fileList[0],unpack=True)
        #print data, data.shape
        #agn = (data[3] > 30 ) & (data[3] < 40 ) & (data[4] > 30 ) & (data[4] < 40 )
        #data_all = data.T[agn]
        #print data_all.shape
        for path_2_input in fileList:
                print path_2_input
                path_2_output = os.path.join(os.environ[env], "light-cone", os.path.basename(path_2_input)[:-4]+".erosita-agn-window-100deg2.gawk.ascii")
                gawk_command = """gawk ' {if ( $4 >= -5 && $4<=5 && $5 >= -5 && $5 <=5 ) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' """ + path_2_input +" > " + path_2_output 
                print(gawk_command)
                os.system(gawk_command)
                #data = n.loadtxt(fileName,unpack=True)
                #print data.shape
                ## compute luminosity
                #dL_cm = (cosmoMD.luminosity_distance(data[2]).to(u.cm)).value
                #flux = 10**(data[9]-0.3) / (4.* n.pi * dL_cm * dL_cm)
                #print dL_cm, flux
                #agn = (flux > 1e-15 ) #& (data[2] < 2.44)
                #print len(agn.nonzero()[0])
                ##data_all = n.vstack((data_all, data.T[agn]))
                ##print data_all.shape
                #n.savetxt(, data.T[agn])    

get_AGN_catalog(env='MD10')
print time.time()-t0, "seconds"

os.system("""cat header_agn.txt MDPL2_ROCKSTAR_FluxProj_*p_000_AGN.erosita-agn-window-100deg2.gawk.ascii > AGN.erosita-agn-window-100deg2.ascii""")
#os.system("""cat  AGN.erosita-agn-window-100deg2.gawk.ascii >  AGN.erosita-agn-window-100deg2-withHeader.ascii""")
