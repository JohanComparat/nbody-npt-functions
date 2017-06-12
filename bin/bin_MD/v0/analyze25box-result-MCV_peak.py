from MultiDark import *
box = MultiDark(Lbox=2500.0 * uu.Mpc,wdir="/data2/DATA/eBOSS/Multidark-lightcones", boxDir = "MD_2.5Gpc",snl =  n.array(glob.glob(join("/data2/DATA/eBOSS/Multidark-lightcones" , "MD_2.5Gpc" , "hlist_?.?????.list"))) ,zsl = None,zArray = n.arange(0.2,2.4,1e-1),Hbox = 67.77 * uu.km / (uu.s * uu.Mpc),Melement = 23593750000.0)


massB=n.arange(8,16,0.01)
vcirB=n.arange(0,4.5,0.01)
concB=n.arange(-5e4,5e4+1,1e3)
n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/Mpeak.bins",massB)
n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/Vpeak.bins",vcirB)
n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/AccRate.bins",concB)

for ii in range(len(box.get_snl())):
    print box.get_snl()[ii].split('_')[-1][:-5]

    centralList = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/*"+box.get_snl()[ii].split('_')[-1][:-5]+"*MVAmatrixCentral*.pkl"))
    centralList.sort()

    nnM = n.empty( [len(centralList),len(massB)-1] ) 
    nnV = n.empty( [len(centralList),len(vcirB)-1] )
    dataVC = n.empty( [len(centralList),len(vcirB)-1,len(concB)-1] )
    dataMC = n.empty( [len(centralList),len(massB)-1,len(concB)-1] )

    for jj in range(len(centralList)):
        f=open(centralList[jj],'r')
        nnMinter,nnVinter,nnCinter,dataMCinter,dataVCinter = cPickle.load(f)
        nnM[jj] = nnMinter
        nnV[jj] = nnVinter 
        dataMC[jj] = dataMCinter[0]
        dataVC[jj] = dataVCinter[0]
        f.close()


    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-central-Mpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([massB[:-1], massB[1:],  nnM.sum(axis=0)]))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-central-Vpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([vcirB[:-1], vcirB[1:],  nnV.sum(axis=0)]) )

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-central-Mpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataMC.sum(axis=0))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-central-Vpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataVC.sum(axis=0))


    satList = n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/*"+box.get_snl()[ii].split('_')[-1][:-5]+"*MVAmatrixSatellite*.pkl"))
    satList.sort()

    nnM = n.empty( [len(satList),len(massB)-1] ) 
    nnV = n.empty( [len(satList),len(vcirB)-1] )
    dataVC = n.empty( [len(satList),len(vcirB)-1,len(concB)-1] )
    dataMC = n.empty( [len(satList),len(massB)-1,len(concB)-1] )

    for jj in range(len(satList)):
        f=open(satList[jj],'r')
        nnMinter,nnVinter,nnCinter,dataMCinter,dataVCinter = cPickle.load(f)
        nnM[jj] = nnMinter
        nnV[jj] = nnVinter 
        dataMC[jj] = dataMCinter[0]
        dataVC[jj] = dataVCinter[0]
        f.close()


    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-sat-Mpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([massB[:-1], massB[1:],  nnM.sum(axis=0)]))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist-sat-Vpeak-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",n.transpose([vcirB[:-1], vcirB[1:],  nnV.sum(axis=0)]) )

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-sat-Mpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataMC.sum(axis=0))

    n.savetxt("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/properties/hist2d-sat-Vpeak-Accrate-"+box.get_snl()[ii].split('_')[-1][:-5]+".dat",dataVC.sum(axis=0))



