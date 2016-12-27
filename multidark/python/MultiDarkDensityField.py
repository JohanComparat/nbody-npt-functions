
"""
.. class:: MultiDark

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class MultiDark is a wrapper to handle Multidark simulations results / outputs.

"""
import cPickle
import fileinput
import astropy.cosmology as co
import astropy.units as u
c2 = co.Planck13
from scipy.interpolate import interp1d
from os.path import join
import os
import astropy.units as uu
import numpy as n
import glob

class MultiDarkSimulation :
    """
    Loads the environement proper to the Multidark simulations. This is the fixed framework of the simulation.
            
    :param Lbox: length of the box in Mpc/h 
    :param wdir: Path to the multidark lightcone directory
    :param boxDir: box directory name
    :param snl: list of snapshots available
    :param zsl: list of redshift corresponding to the snapshots   
    :param zArray: redshift array to be considered to interpolate the redshift -- distance conversion
    :param Hbox: Hubble constant at redshift 0 of the box
    :param Melement: Mass of the resolution element in solar masses.   
    :param columnDict: dictionnary to convert column name into the index to find it in the snapshots
    """

    def __init__(self,Lbox=2500.0 * uu.Mpc, wdir="/data2/DATA/eBOSS/Multidark-lightcones/", boxDir="MD_2.5Gpc", snl=n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_?.?????.list")), zsl=None, zArray=n.arange(0.2,2.4,1e-1), Hbox = 67.77 * uu.km / (uu.s * uu.Mpc), Melement = 23593750000.0 ):
        self.Lbox = Lbox # box length
        self.Hbox = Hbox # Hubble constant at redshift 0 in the box
        self.wdir = wdir # working directory
        self.boxDir = boxDir # directory of the box where the snapshots a stored
        self.snl = snl # snapshot list
        self.zsl = zsl # corresponding redshift list
        self.zArray = zArray # redshift for the dC - z conversion
        self.Melement = Melement # mass of one particle in the box
        self.h = 0.6777
        self.omega_lambde = 0.692885
        self.omega_matter = 0.307115
        self.omega_baryon = 0.048206
        self.ns = 0.96
        self.sigma8 = 0.8228
        self.G = 6.67428 * 10**(-9) # cm3 g-1 s-2
        self.Msun = 1.98892 * 10**(33.) # g
        self.Npart = 3840
        self.force_resolution = 5. # kpc /h
        self.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np' : 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin': 17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23,'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'c_to_a_500c': 33, 'Ax_500c': 34, 'Ay_500c': 35, 'Az_500c': 36, 'TU': 37, 'M_pe_Behroozi': 38, 'M_pe_Diemer': 39, 'pid': 40}

        if self.boxDir == "MD_0.4Gpc":
            self.Melement = 9.63 * 10**7 # Msun
            self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')

        if self.boxDir == "MD_1Gpc_new_rockS":
            self.Melement = 1.51 * 10**9. # Msun
            self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')

        if self.boxDir == "MD_2.5Gpc":
            self.Melement = 2.359 * 10**10. # Msun
            self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')

        if self.boxDir == "MD_4Gpc":
            self.Melement = 9.6 * 10**10. # Msun
            self.Npart = 4096
            self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')

        if self.boxDir == "MD_2.5Gpc":
            # for satellites ...
            self.columnDictHlist = {'scale': 0, 'id': 1, 'desc_scale': 2, 'desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Rs_Klypin': 35, 'Mmvir_all': 36, 'M200b': 37, 'M200c': 38, 'M500c': 39, 'M2500c': 40, 'Xoff': 41, 'Voff': 42, 'Spin_Bullock': 43, 'b_to_a': 44, 'c_to_a': 45, 'Ax': 46, 'Ay': 47, 'Az': 48, 'b_to_a_500c': 49, 'c_to_a_500c': 50, 'Ax_500c': 51, 'Ay_500c': 52, 'Az_500c': 53, 'TU': 54, 'M_pe_Behroozi': 55, 'M_pe_Diemer': 56, 'Halfmass_Radius': 57, 'Macc': 58, 'Mpeak': 59, 'Vacc': 60, 'Vpeak': 61, 'Halfmass_Scale': 62, 'Acc_Rate_Inst': 63, 'Acc_Rate_100Myr': 64, 'Acc_Rate_1Tdyn': 65, 'Acc_Rate_2Tdyn': 66, 'Acc_Rate_Mpeak': 67, 'Mpeak_Scale': 68, 'Acc_Scale': 69, 'First_Acc_Scale': 70, 'First_Acc_Mvir': 71, 'First_Acc_Vmax': 72, 'VmaxatMpeak': 73}
        if self.boxDir == "MD_0.4Gpc" or self.boxDir == "MD_1Gpc_new_rockS" :
            # for satellites ...
            self.columnDictHlist = {'scale': 0, 'id': 1, 'desc_scale': 2, 'desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Rs_Klypin': 35, 'Mmvir_all': 36, 'M200b': 37, 'M200c': 38, 'M500c': 39, 'M2500c': 40, 'Xoff': 41, 'Voff': 42, 'Spin_Bullock': 43, 'b_to_a': 44, 'c_to_a': 45, 'Ax': 46, 'Ay': 47, 'Az': 48, 'b_to_a_500c': 49, 'c_to_a_500c': 50, 'Ax_500c': 51, 'Ay_500c': 52, 'Az_500c': 53, 'TU': 54, 'M_pe_Behroozi': 55, 'M_pe_Diemer': 56, 'Macc': 57, 'Mpeak': 58, 'Vacc': 59, 'Vpeak': 60, 'Halfmass_Scale': 61, 'Acc_Rate_Inst': 62, 'Acc_Rate_100Myr': 63, 'Acc_Rate_1Tdyn': 64, 'Acc_Rate_2Tdyn': 65, 'Acc_Rate_Mpeak': 66, 'Mpeak_Scale': 67, 'Acc_Scale': 68, 'First_Acc_Scale': 69, 'First_Acc_Mvir': 70, 'First_Acc_Vmax': 71, 'VmaxatMpeak': 72}

	def get_DF_at_XYZ(x, y, z, path_to_DF, Lbox=1000., gridSize = 2048.):
		dL = Lbox/gridSize
		sel =( x > dL*(0.5 + ii) ) & ( x < dL*(0.5 + ii + 1) ) & ( y > dL*(0.5 + jj) ) & ( y < dL*(0.5 + jj + 1) ) & ( x > dL*(0.5 + kk) ) & ( z < dL*(0.5 + kk + 1) )
		compute : 
		imax = x/dL - 0.5 
		imin = x/dL - 0.5 - 1 
		jmax = y/dL - 0.5 
		jmin = y/dL - 0.5 - 1 
		kmax = z/dL - 0.5 
		kmin = z/dL - 0.5 - 1 
		f=open(path_to_DF,'r')
		qty = n.empty( (Nratio,len(bins)-1) ) 
		data1 =  n.fromfile(f,dtype="float64",count=NperBatch) # 512 cube

		
    def computeSingleDistributionFunction(self, ii, name, bins, Mfactor=100. ) :
        """
        Extracts the distribution of quantity 'name' out of all snapshots of the Multidark simulation.        
        :param ii: index of the snapshot in the list self.snl
        :param name: name of the quantity of interest, mass, velocity.
        :param index: of the quantity of interest in the snapshots.
        :param bins: binning scheme to compute the historgram.
        :param Mfactor: only halos with Mvir > Mfact* Melement are used.
        """		
        index = self.columnDict[name]
        output_dir = join(self.wdir,self.boxDir,"properties",name)
        os.system('mkdir '+ output_dir)
        NperBatch = 10000000
        qtyCentral = n.empty(NperBatch)  # 10M array
        qtySat = n.empty(NperBatch)  # 10M array
        print name, index, output_dir

        fl = fileinput.input(self.snl[ii])
        nameSnapshot = self.snl[ii].split('/')[-1][:-5]

        countCen,countSat,countFileCen,countFileSat = 0,0,0,0
        
        for line in fl:
            if line[0] == "#" :
                continue

            line = line.split()
            sat_or_cen = float(line[self.columnDict['pid']])
            mv = float(line[self.columnDict['mvir']])
            if sat_or_cen != -1 and mv > Mfactor * self.Melement :
                countSat+= 1					
                qtySat[countSat] = float(line[index])
                
            if sat_or_cen == -1 and mv > Mfactor * self.Melement :
                countCen+= 1					
                qtyCentral[countCen] = float(line[index])
                
            if countCen == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(qtyCentral),bins = bins)
                print "countCen",countCen
                f = open(join(output_dir, nameSnapshot + "_" + name + "_Central_" + str(countFileCen)+ ".pkl"),'w')
                cPickle.dump(nnM,f)
                f.close()
                countFileCen+= 1
                countCen = 0
                qtyCentral = n.empty(NperBatch)

            if countSat == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(qtySat),bins = bins)
                print "countSat", countSat
                f = open(join(output_dir, nameSnapshot + "_" + name+ "_Satellite_" + str(countFileSat)+ ".pkl"),'w')
                cPickle.dump(nnM,f)
                f.close()
                countFileSat+= 1
                countSat = 0
                qtySat = n.empty(NperBatch)

        # and for the last batch :
        nnM,bb = n.histogram(n.log10(qtyCentral),bins = bins)
        f = open(join(output_dir, nameSnapshot + "_" + name +"_Central_" + str(countFileCen)+ ".pkl"),'w')
        cPickle.dump(nnM,f)
        f.close()

        nnM,bb = n.histogram(n.log10(qtySat),bins = bins)
        f = open(join(output_dir, nameSnapshot + "_" + name + "_Satellite_" + str(countFileSat)+ ".pkl"),'w')
        cPickle.dump(nnM,f)
        f.close()
        
        n.savetxt(join(output_dir,name+".bins"),n.transpose([bins]))


    def computeDensityFieldDistributionFunction(self, path_to_DF, outputFile, bins ) :
		"""
		Extracts the distribution of quantity 'name' out of all snapshots of the Multidark simulation.        
		:param path_to_DF: path to the density field file
		:param outputFile: where the histogram iswritten
		:param bins: binning scheme to compute the historgram.
		"""		
		NperBatch = 512**3
		Ntotal = 2048**3
		Nratio = Ntotal / NperBatch
		#qtyCentral = n.empty( (64,NperBatch) ) 
		f=open(path_to_DF,'r')
		qty = n.empty( (Nratio,len(bins)-1) ) 
		for ii in n.arange(Nratio):
			data1 =  n.fromfile(f,dtype="float64",count=NperBatch) # 512 cube
			nnM,bb = n.histogram(n.log10(data1),bins = bins)
			qty[ii] = nnM
				
		n.savetxt(outputFile+".hist",n.transpose([bins[:-1], bins[1:], qty.sum(axis=0) ]), header = " minLogDelta maxLogDelta N")

		
    def computeDensityFieldForHaloNumber(self, path_to_DF, path_to_RS, outputFile, gridSize=2048, subgridSize = 256 ) :
		"""
        Extracts the distribution of quantity 'name' out of all snapshots of the Multidark simulation.        
        :param path_to_DF: path to the density field file
        :param path_to_RS: path to the rockstar halo catalog file
        :param outputFile: where the histogram iswritten
        :param bins: binning scheme to compute the historgram.
		:param gridSize: grid size from the density field
		:param subgridSize: grid size to compute histograms on and write outputs.
        """		
        #In each cell, average the number of counts in the RS file : Ncen Nsat
		dL = 1000./gridSize
		NperBatch = subgridSize**3
		Ntotal = gridSize**3
		Nratio = Ntotal / NperBatch

		hf=open(path_to_RS,'r')
		DEFINE x, y, z, c_o_s, 
		REWRITE halos of interest ?
		
		f=open(path_to_DF,'r')
		out = n.empty( (subgridSize**3, 3) )
		count = 0
		countOut = 0
		for kk in n.arange(gridSize):
			for jj in n.arange(gridSize):
				for ii in n.arange(gridSize):
					sel =( x > dL*(0.5 + ii) ) & ( x < dL*(0.5 + ii + 1) ) & ( y > dL*(0.5 + jj) ) & ( y < dL*(0.5 + jj + 1) ) & ( x > dL*(0.5 + kz) ) & ( z < dL*(0.5 + kk + 1) )
					Nhalos = len(sel.nonzero()[0])
					selCen =  (sel) & (CONDITION_CEN)
					NCentrals = len(selCen.nonzero()[0])
					deltaValue = n.fromfile(f,dtype="float64",count=1)
					out[count] = n.array([deltaValue, Nhalos, Ncentrals])
					if count == subgridSize**3 :
						dataAB_tot = n.histogram2d(n.log10(out.T[0]), n.log10(out.T[1]) ,bins = [binsA,binsB])
						dataAB_cen = n.histogram2d(n.log10(out.T[0]), n.log10(out.T[2]) ,bins = [binsA,binsB])
						f = open(outputFile + "_" +str(countOut)+".pkl" ,'w')
						cPickle.dump([binsA,binsB,dataAB_tot, dataAB_cen],f)
						f.close()

						out = n.empty( (subgridSize**3, 3) )
						countOut +=1 
					
					count += 1
					
		#GATHER RESULTS
		fileList = glob.glob(outputFile + "_*.pkl")
		out_all = n.empty( (Nratio, len(binsA)-1, len(binsB)-1) )
		out_cen = n.empty( (Nratio, len(binsA)-1, len(binsB)-1) )
		for ii, el in enumerate(fileList):
			f = open(el ,'r')
			binsA,binsB,dataAB_tot, dataAB_cen = cPickle.dump([binsA,binsB,dataAB_tot, dataAB_cen],f)
			f.close()
			out_all[ii] = dataAB_tot
			out_cen[ii] = dataAB_cen

		
		f = open(outputFile + "_all.pkl" ,'w')
		cPickle.dump([binsA,binsB,n.sum(out_all, xis=0), n.sum(out_cen, xis=0)],f)
		f.close()

		
    def computeDensityFieldHaloCorrelation(self, path_to_DF, path_to_RS, outputFile, bins ) :
		"""
        Extracts the distribution of quantity 'name' out of all snapshots of the Multidark simulation.        
        :param path_to_DF: path to the density field file
        :param path_to_RS: path to the rockstar halo catalog file
        :param outputFile: where the histogram iswritten
        :param bins: binning scheme to compute the historgram.
        """		
        #In each cell, average the number of counts in the RS file : Ncen Nsat
		NperBatch = 512**3
		Ntotal = 2048**3
		Nratio = Ntotal / NperBatch

		dL = 1000/2048.
		
		f=open(path_to_DF,'r')
		for N in n.arange(Ntotal):
			deltaValue = n.fromfile(f,dtype="float64",count=1) # 512 cube
			sel =( x > dL*(0.5 + N) ) & ( x < dL*(0.5 + N + 1) )
		

    def combinesSingleDistributionFunction(self, ii, name='Vpeak', bins=10**n.arange(0,3.5,0.01), type = "Central" ) :
        """
        Coombines the outputs of computeSingleDistributionFunction.
        :param ii: index of the snapshot
        :param name: name of the quantity studies
        :param bins: bins the histogram was done with
        :param type: "Central" or "Satellite"
        """
        output_dir = join(self.wdir,self.boxDir,"properties",name)
        nameSnapshot = self.snl[ii].split('/')[-1][:-5]
        pklList = n.array(glob.glob(join(output_dir, nameSnapshot + "_" + name +"_"+type+"_*.pkl")))

        nnM = n.empty( [len(pklList),len(bins)-1] ) 
        for jj in range(len(pklList)):
            f=open(pklList[jj],'r')
            nnMinter = cPickle.load(f)
            nnM[jj] = nnMinter
            f.close()

        n.savetxt(join(output_dir,"hist-"+type+"-"+name+"-"+nameSnapshot[6:]+".dat"),n.transpose([bins[:-1], bins[1:], nnM.sum(axis=0)]))


    def computeDoubleDistributionFunction(self, ii, nameA, nameB, binsA, binsB, Mfactor = 100. ) :
        """
        Extracts the distributions of two quantity and their correlation 'name' out of all snapshots of the Multidark simulation.
        :param ii: index of the snapshot in the list self.snl
        :param name: name of the quantity of interest, mass, velocity.
        :param index: of the quantity of interest in the snapshots.
        :param bins: binning scheme to compute the historgram.
        :param Mfactor: only halos with Mvir > Mfact* Melement are used.
        """		
        indexA = self.columnDict[nameA]
        indexB = self.columnDict[nameB]
        output_dir = join(self.wdir,self.boxDir,"properties",nameA+"-"+nameB)
        os.system('mkdir '+ output_dir)
        NperBatch = 10000000
        qtyCentral = n.empty((NperBatch,2))  # 10M array
        qtySat = n.empty((NperBatch,2))  # 10M array
        print nameA, nameB, indexA, indexB, output_dir

        fl = fileinput.input(self.snl[ii])
        nameSnapshot = self.snl[ii].split('/')[-1][:-5]

        countCen,countSat,countFileCen,countFileSat = 0,0,0,0
        
        for line in fl:
            if line[0] == "#" :
                continue

            line = line.split()
            sat_or_cen = float(line[self.columnDict['pid']])
            mv = float(line[self.columnDict['mvir']])
            if sat_or_cen != -1 and mv > Mfactor * self.Melement :
                countSat+= 1					
                qtySat[countSat] = float(line[indexA]),float(line[indexB])
                
            if sat_or_cen == -1 and mv > Mfactor * self.Melement :
                countCen+= 1					
                qtyCentral[countCen] = float(line[indexA]),float(line[indexB])
                
            if countCen == NperBatch-1 :
                nnA,bbA = n.histogram(n.log10(qtyCentral.T[0]),bins = binsA)
                nnB,bbB = n.histogram(n.log10(qtyCentral.T[1]),bins = binsB)
                dataAB = n.histogram2d(n.log10(qtyCentral.T[0]), n.log10(qtyCentral.T[1]) ,bins = [binsA,binsB])
                print "countCen",countCen
                f = open(join(output_dir, nameSnapshot + "_" + nameA+"-"+nameB + "_Central_" + str(countFileCen)+ ".pkl"),'w')
                cPickle.dump([nnA,nnB,dataAB],f)
                f.close()
                countFileCen+= 1
                countCen = 0
                qtyCentral = n.empty((NperBatch,2))

            if countSat == NperBatch-1 :
                nnA,bbA = n.histogram(n.log10(qtySat.T[0]),bins = binsA)
                nnB,bbB = n.histogram(n.log10(qtySat.T[1]),bins = binsB)
                dataAB = n.histogram2d(n.log10(qtySat.T[0]), n.log10(qtySat.T[1]) ,bins = [binsA,binsB])
                print "countSat", countSat
                f = open(join(output_dir, nameSnapshot + "_" + nameA+"-"+nameB+ "_Satellite_" + str(countFileSat)+ ".pkl"),'w')
                cPickle.dump([nnA,nnB,dataAB],f)
                f.close()
                countFileSat+= 1
                countSat = 0
                qtySat = n.empty((NperBatch,2))

        # and for the last batch :
        nnA,bbA = n.histogram(n.log10(qtyCentral.T[0]),bins = binsA)
        nnB,bbB = n.histogram(n.log10(qtyCentral.T[1]),bins = binsB)
        dataAB = n.histogram2d(n.log10(qtyCentral.T[0]), n.log10(qtyCentral.T[1]) ,bins = [binsA,binsB])
        print "countCen",countCen
        f = open(join(output_dir, nameSnapshot + "_" + nameA+"-"+nameB + "_Central_" + str(countFileCen)+ ".pkl"),'w')
        cPickle.dump([nnA,nnB,dataAB],f)
        f.close()

        nnA,bbA = n.histogram(n.log10(qtySat.T[0]),bins = binsA)
        nnB,bbB = n.histogram(n.log10(qtySat.T[1]),bins = binsB)
        dataAB = n.histogram2d(n.log10(qtySat.T[0]), n.log10(qtySat.T[1]) ,bins = [binsA,binsB])
        print "countSat", countSat
        f = open(join(output_dir, nameSnapshot + "_" + nameA+"-"+nameB+ "_Satellite_" + str(countFileSat)+ ".pkl"),'w')
        cPickle.dump([nnA,nnB,dataAB],f)
        f.close()
        
        n.savetxt(join(output_dir,nameA+".bins"),n.transpose([binsA]))
        n.savetxt(join(output_dir,nameB+".bins"),n.transpose([binsB]))



    def combinesDoubleDistributionFunction(self, ii, nameA, nameB, binsA, binsB, type = "Central" ) :
        """
        Coombines the outputs of computeDoubleDistributionFunction.
        :param ii: index of the snapshot
        :param name: name of the quantity studies
        :param bins: bins the histogram was done with
        :param type: "Central" or "Satellite"
        """
        output_dir = join(self.wdir,self.boxDir,"properties",nameA+"-"+nameB)
        nameSnapshot = self.snl[ii].split('/')[-1][:-5]
        pklList = n.array(glob.glob(join(output_dir, nameSnapshot + "_" + nameA+"-"+nameB +"_"+type+"_*.pkl")))

        nnA = n.empty( [len(pklList),len(binsA)-1] ) 
        nnB = n.empty( [len(pklList),len(binsB)-1] ) 
        dataAB = n.empty( [len(pklList),len(binsA)-1,len(binsB)-1] ) 
        for jj in range(len(pklList)):
            f=open(pklList[jj],'r')
            nnAinter, nnBinter, dataABinter = cPickle.load(f)
            nnA[jj] = nnAinter
            nnB[jj] = nnBinter
            dataAB[jj] = dataABinter[0]
            f.close()

        n.savetxt(join(output_dir,"hist-"+type+"-"+nameA+"-"+nameSnapshot[6:]+".dat"),n.transpose([binsA[:-1], binsA[1:], nnA.sum(axis=0)]))
        n.savetxt(join(output_dir,"hist-"+type+"-"+nameB+"-"+nameSnapshot[6:]+".dat"),n.transpose([binsB[:-1], binsB[1:], nnB.sum(axis=0)]))
        n.savetxt(join(output_dir, "hist2d-"+type+"-"+ nameA+"-"+nameB + "-"+ nameSnapshot[6:] + ".dat"), dataAB.sum(axis=0))


    def computeMassVelocityConcentrationFunction(self,ii) :
        """
        DO NOT USE
        computes the mass, velocity and concentration histograms for a rockstar snapshot.
        :param ii: index of the snapshot in the list self.snl
        # does not work any more 
        DO NOT USE
        """
        massB = n.arange(8,16,0.01)
        vcirB = n.arange(0,4.5,0.01)
        concB = n.arange(1,3,0.1)

        NperBatch = 10000000
        mvcCentralMatrix = n.empty((NperBatch,3))  # 1M matrixes
        mvcSatMatrix = n.empty((NperBatch,3))  # 1 M matrixes

        fl = fileinput.input(self.snl[ii])
        name = self.snl[ii].split('/')[-1][:-5]
        countCen,countSat,countFileCen,countFileSat = 0,0,0,0
        for line in fl:
            if line[0] == "#" :
                continue

            line = line.split()
            sat_or_cen = float(line[5])
            if sat_or_cen != -1 :
                countSat+= 1					
                mvcSatMatrix[countSat] = float(line[10]), float(line[16]), float(line[11]) 
                
            if sat_or_cen == -1 :
                countCen+= 1					
                mvcCentralMatrix[countCen] = float(line[10]), float(line[16]), float(line[11])
                
            if countCen == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
                nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
                nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
                dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
                dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
                print "countCen",countCen
                f = open(join(self.wdir,self.boxDir,"properties", name+"_MVRmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
                cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
                f.close()
                countFileCen+= 1
                countCen = 0

            if countSat == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
                nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
                nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
                dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
                dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
                print "countSat", countSat
                f = open(join(self.wdir,self.boxDir ,"properties" , 
    name+"_MVRmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
                cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
                f.close()
                countFileSat+= 1
                countSat = 0

        # and for the last batch :
        nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
        nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
        nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
        dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
        dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
        f = open(join(self.wdir,self.boxDir,"properties",name+ "_MVRmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
        cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
        f.close()

        nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
        nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
        nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
        dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
        dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
        f = open(join(self.wdir,self.boxDir,"properties",name+ "_MVRmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
        cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
        f.close()


    def computeMassVelocityPeakAccRateFunctions(self,ii) :
        """
        DO NOT USE
        computes the mass, velocity and concentration histograms for a rockstar snapshot.
        :param ii: index of the snapshot in the list self.snl()
        DO NOT USE
        """
        massB = n.arange(8,16,0.01)
        vcirB = n.arange(0,4.5,0.01)
        concB = n.arange(-5e4,5e4+1,1e3)

        NperBatch = 10000000
        mvcCentralMatrix = n.empty((NperBatch,3))  # 1M matrixes
        mvcSatMatrix = n.empty((NperBatch,3))  # 1 M matrixes

        fl = fileinput.input(self.snl[ii])
        name = self.snl[ii].split('/')[-1][:-5]
        countCen,countSat,countFileCen,countFileSat = 0,0,0,0
        for line in fl:
            if line[0] == "#" :
                continue

            line = line.split()
            sat_or_cen = float(line[5])
            if sat_or_cen != -1 :
                countSat+= 1					
                #print mvcSatMatrix[countSat]
                #print line[59], line[61], line[67]
                mvcSatMatrix[countSat] = float(line[59]), float(line[61]), float(line[67]) # check the right indices ... MASS velocity concentration

            if sat_or_cen == -1 :
                countCen+= 1					
                #print mvcCentralMatrix[countCen]
                #print line[59], line[61], line[67]
                mvcCentralMatrix[countCen] = float(line[59]), float(line[61]), float(line[67]) # check the right indices ... MASS velocity concentration

            if countCen == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
                nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
                nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
                dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
                dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
                print "countCen",countCen
                f = open(join(self.wdir,self.boxDir,"properties", name+"_MVAmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
                cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
                f.close()
                countFileCen+= 1
                countCen = 0

            if countSat == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
                nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
                nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
                dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
                dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
                print "countSat", countSat
                f = open(join(self.wdir,self.boxDir ,"properties" , 
    name+"_MVAmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
                cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
                f.close()
                countFileSat+= 1
                countSat = 0

        # and for the last batch :
        nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
        nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
        nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
        dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
        dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
        f = open(join(self.wdir,self.boxDir,"properties",name+ "_MVAmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
        cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
        f.close()

        nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
        nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
        nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
        dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
        dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
        f = open(join(self.wdir,self.boxDir,"properties",name+ "_MVAmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
        cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
        f.close()


