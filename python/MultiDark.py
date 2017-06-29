
"""
.. class:: MultiDark

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class MultiDark is a wrapper to handle Multidark simulations results / outputs.

"""
from os.path import join
import os
import glob
import time

import cPickle
import fileinput
import astropy.io.fits as fits

import numpy as n
from scipy.interpolate import interp1d
import scipy.spatial.ckdtree as t

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115, Ob0=0.048206)
cosmoDS = FlatLambdaCDM(H0=68.46*u.km/u.s/u.Mpc, Om0=0.298734, Ob0=0.046961)
import astropy.constants as constants

G =  constants.G.to(u.kpc**3/(u.solMass * u.yr**2)).value

t_dynamical = lambda rvir, mvir : (rvir**3./(G*mvir))**0.5
def tau_quenching(tdyn, tau_0, tau_s, m_star):
	if m_star < 1e10 :
		return tdyn * tau_0
	else :
		return tdyn * tau_0 * (m_star * 10.**(-10.))**(tau_s)

f_loss = lambda t : 0.05*n.ln( 1 + t / (1.4*10**6))


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

	def __init__(self,Lbox=2500.0 * u.Mpc, wdir="/data2/DATA/eBOSS/Multidark-lightcones/", boxDir="MD_2.5Gpc", snl=n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_?.?????.list")), zsl=None, zArray=n.arange(0.2,2.4,1e-1), Hbox = 67.77 * u.km / (u.s * u.Mpc), Melement = 23593750000.0 ):
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
		self.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'c_to_a_500c': 33, 'Ax_500c': 34, 'Ay_500c': 35, 'Az_500c': 36, 'TU': 37, 'M_pe_Behroozi': 38, 'M_pe_Diemer': 39, 'pid': 40}
		#self.columnDictHlist = {'scale': 0, 'id': 1, 'desc_scale': 2, 'desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Rs_Klypin': 35, 'Mmvir_all': 36, 'M200b': 37, 'M200c': 38, 'M500c': 39, 'M2500c': 40, 'Xoff': 41, 'Voff': 42, 'Spin_Bullock': 43, 'b_to_a': 44, 'c_to_a': 45, 'Ax': 46, 'Ay': 47, 'Az': 48, 'b_to_a500c': 49, 'c_to_a500c': 50, 'Ax500c': 51, 'Ay500c': 52, 'Az500c': 53, 'TU': 54, 'M_pe_Behroozi': 55, 'M_pe_Diemer': 56, 'Macc': 57, 'Mpeak': 58, 'Vacc': 59, 'Vpeak': 60, 'Halfmass_Scale': 61, 'Acc_Rate_Inst': 62, 'Acc_Rate_100Myr': 63, 'Acc_Rate_1Tdyn': 64, 'Acc_Rate_2Tdyn': 65, 'Acc_Rate_Mpeak': 66, 'Mpeak_Scale': 67, 'Acc_Scale': 68, 'First_Acc_Scale': 69, 'First_Acc_Mvir': 70, 'First_Acc_Vmax': 71, 'VmaxAtMpeak': 72}
		self.columnDictHlist = {'scale': 0, 'id': 1, 'desc_scale': 2, 'desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Tidal_Force': 35, 'Tidal_ID': 36, 'Rs_Klypin': 37, 'Mmvir_all': 38, 'M200b': 39, 'M200c': 40, 'M500c': 41, 'M2500c': 42, 'Xoff': 43, 'Voff': 44, 'Spin_Bullock': 45, 'b_to_a': 46, 'c_to_a': 47, 'Ax': 48, 'Ay': 49, 'Az': 50, 'b_to_a_500c' : 51, 'c_to_a_500c' : 52, 'Ax_500c' : 53, 'Ay_500c' : 54, 'Az_500c' : 55, 'TU': 56, 'M_pe_Behroozi': 57, 'M_pe_Diemer': 58, 'Macc': 59, 'Mpeak': 60, 'Vacc': 61, 'Vpeak': 62, 'Halfmass_Scale': 63, 'Acc_Rate_Inst': 64, 'Acc_Rate_100Myr': 65, 'Acc_Rate_1Tdyn': 66, 'Acc_Rate_2Tdyn': 67, 'Acc_Rate_Mpeak': 68, 'Mpeak_Scale': 69, 'Acc_Scale': 70, 'First_Acc_Scale': 71, 'First_Acc_Mvir': 72, 'First_Acc_Vmax': 73, 'VmaxAtMpeak': 74, 'Tidal_Force_Tdyn': 75, 'logVmaxVmaxmaxTdynTmpeak': 76, 'Time_to_future_merger': 77, 'Future_merger_MMP_ID': 78 }
		
		self.columnDictHlist25 = {'scale': 0, 'id': 1, 'desc_scale': 2, 'desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Rs_Klypin': 35, 'Mmvir_all': 36, 'M200b': 37, 'M200c': 38, 'M500c': 39, 'M2500c': 40, 'Xoff': 41, 'Voff': 42, 'Spin_Bullock': 43, 'b_to_a': 44, 'c_to_a': 45, 'Ax': 46, 'Ay': 47, 'Az': 48, 'b_to_a_500c' : 49, 'c_to_a_500c' : 50, 'Ax_500c' : 51, 'Ay_500c' : 52, 'Az_500c' : 53, 'TU': 54, 'M_pe_Behroozi': 55, 'M_pe_Diemer': 56, 'Halfmass_Radius': 57, 'Macc': 58, 'Mpeak': 59, 'Vacc': 60, 'Vpeak': 61, 'Halfmass_Scale': 62, 'Acc_Rate_Inst': 63, 'Acc_Rate_100Myr': 64, 'Acc_Rate_1Tdyn': 65, 'Acc_Rate_2Tdyn': 66, 'Acc_Rate_Mpeak': 67, 'Mpeak_Scale': 68, 'Acc_Scale': 69, 'First_Acc_Scale': 70, 'First_Acc_Mvir': 71, 'First_Acc_Vmax': 72, 'VmaxAtMpeak': 73, 'Tidal_Force_Tdyn': 74, 'logVmaxVmaxmaxTdynTmpeak': 75, 'Time_to_future_merger': 76, 'Future_merger_MMP_ID': 77 }
		
		if self.boxDir == "MD_0.4Gpc":
			self.Melement = 9.63 * 10**7 # Msun
			#self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')
			
		if self.boxDir == "MD_0.4GpcNW":
			self.Melement = 9.63 * 10**7 # Msun
			#self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')

		if self.boxDir == "MD_1Gpc" or self.boxDir == "MD_1.0Gpc":
			self.Melement = 1.51 * 10**9. # Msun
			#self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')

		if self.boxDir == "MD_2.5Gpc":
			self.Melement = 2.359 * 10**10. # Msun
			#self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')
			
		if self.boxDir == "MD_2.5GpcNW":
			self.Melement = 2.359 * 10**10. # Msun
			self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')
			self.columnDict = {'id': 0, 'desc_id': 1, 'mvir': 2, 'vmax': 3, 'vrms': 4, 'rvir': 5, 'rs': 6, 'Np': 7, 'x': 8, 'y': 9, 'z': 10, 'vx': 11, 'vy': 12, 'vz': 13, 'Jx': 14, 'Jy': 15, 'Jz': 16, 'Spin':17, 'Rs_Klypin': 18, 'Mmvir_all': 19, 'M200b': 20, 'M200c': 21, 'M500c': 22, 'M2500c': 23, 'Xoff': 24, 'Voff': 25, 'Spin_Bullock': 26, 'b_to_a': 27, 'c_to_a': 28, 'Ax': 29, 'Ay': 30, 'Az': 31, 'b_to_a_500c': 32, 'c_to_a_500c': 33, 'Ax_500c': 34, 'Ay_500c': 35, 'Az_500c': 36, 'TU': 37, 'M_pe_Behroozi': 38, 'M_pe_Diemer': 39, 'Halfmass_Radius': 40, 'pid': 41}

		if self.boxDir == "MD_4Gpc" or self.boxDir == "MD_4.0Gpc":
			self.Melement = 9.6 * 10**10. # Msun
			self.Npart = 4096
			#self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')

		if self.boxDir == "MD_4GpcNW":
			self.Melement = 9.6 * 10**10. # Msun
			self.Npart = 4096
			#self.vmin = 4* (self.Melement*self.Msun*self.G/(self.force_resolution*u.kpc.to('cm')))**0.5 * u.cm.to('km')
			
	def writeEMERGEcatalog(self, path_2_snapshot, rho_crit, delta_vir, mmin=10**8, NperBatch = 2000000, file_identifier = "_EMERGE_Nb_"):
		"""
		Extracts the positions and mass out of a snapshot of the Multidark simulation.        
		:param ii: index of the snapshot in the list self.snl
		:param vmin: name of the quantity of interest, mass, velocity.
		:param vmax: of the quantity of interest in the snapshots.
		:param NperBatch: number of line per fits file, default: 1000000
		"""
		# sets the columns definition
		self.columnDict = self.columnDictHlist 
		# opens the hlist file
		fl = fileinput.input(path_2_snapshot)
		nameSnapshot = os.path.basename(path_2_snapshot)[:-5]
		Nb = 0
		count = 0
		output = n.zeros((NperBatch,12))
		outfile_base = os.path.join(os.environ["MD10"],"emerge",nameSnapshot + file_identifier)
		for line in fl:
			#print line
			if line[0] == "#" :
				continue

			line = line.split()
			# compute tdyn, 
			rvir = float(line[self.columnDict['rvir']])
			mvir = float(line[self.columnDict['mvir']])
			tdyn = t_dynamical(rvir, mvir)
			
			# compute density at rvir
			rs = float(line[self.columnDict['rs']])
			conc = rvir / rs
			rho_at_rvir = rho_crit * delta_vir * conc**2. / ((1+conc)*((1+conc)*n.ln(1+conc)-conc))
			
			newline =n.array([ int(line[self.columnDict['id']]), int(line[self.columnDict['pid']]), int(line[self.columnDict['Snap_num']]), rvir, n.log10(mvir), n.log10(float(line[self.columnDict['Mpeak']])), float(line[self.columnDict['Mpeak_Scale']]), float(line[self.columnDict['Acc_Rate_1Tdyn']]), float(line[self.columnDict['Time_to_future_merger']]), float(line[self.columnDict['Future_merger_MMP_ID']]), tdyn, rho_at_rvir])
			
			
			#print newline
			#print newline.shape
			if float(line[self.columnDict['mvir']])>mmin :
				output[count] = newline
				count+=1
				
			if count == NperBatch  :
				out_filename = outfile_base + str(Nb).zfill(3) + ".fits"
				print( out_filename )
				#print "count",count
				#print output
				#print output.shape
				#print output.T[0].shape
				#define the columns
				hdu_cols  = fits.ColDefs([
				fits.Column(name='id',format='I',            array= output.T[0] )
				,fits.Column(name='pid',format='I',         array= output.T[1] ) 
				,fits.Column(name='Snap_num',format='I',         array= output.T[2] ) 
				,fits.Column(name='rvir',format='D',         array= output.T[3] ) 
				,fits.Column(name='mvir',format='D',           array= output.T[4] ) 
				,fits.Column(name='Mpeak',format='D',          array= output.T[5] ) 
				,fits.Column(name='Mpeak_scale',format='D',            array= output.T[6] ) 
				,fits.Column(name='Acc_Rate_1Tdyn',format='D',            array= output.T[7] ) 
				,fits.Column(name='Time_to_future_merger',format='D',            array= output.T[8] ) 
				,fits.Column(name='Future_merger_MMP_ID',format='D',           array= output.T[9] ) 
				,fits.Column(name='tdyn',format='D',           array= output.T[10] ) 
				,fits.Column(name='rho_at_rvir',format='D',           array= output.T[11] ) 
				])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				#define the header
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = nameSnapshot
				prihdr['count'] = count
				prihdr['batchN'] = Nb
				prihdr['author'] = 'JC'
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				out_filename = outfile_base + str(Nb).zfill(3)+".fits"
				print( out_filename )
				os.system("rm "+out_filename)
				thdulist.writeto(out_filename)
				Nb+=1
				count=0
				#resest the output matrix
				output = n.zeros((NperBatch,12))
		
		
		# and for the last batch :		
		hdu_cols  = fits.ColDefs([
		fits.Column(name='id',format='I',            array= output.T[0] )
		,fits.Column(name='pid',format='I',         array= output.T[1] ) 
		,fits.Column(name='Snap_num',format='I',         array= output.T[2] ) 
		,fits.Column(name='rvir',format='D',         array= output.T[3] ) 
		,fits.Column(name='mvir',format='D',           array= output.T[4] ) 
		,fits.Column(name='Mpeak',format='D',          array= output.T[5] ) 
		,fits.Column(name='Mpeak_scale',format='D',            array= output.T[6] ) 
		,fits.Column(name='Acc_Rate_1Tdyn',format='D',            array= output.T[7] ) 
		,fits.Column(name='Time_to_future_merger',format='D',            array= output.T[8] ) 
		,fits.Column(name='Future_merger_MMP_ID',format='D',           array= output.T[9] ) 
		,fits.Column(name='tdyn',format='D',           array= output.T[10] ) 
		,fits.Column(name='rho_at_rvir',format='D',           array= output.T[11] ) 
		])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		#define the header
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = nameSnapshot
		prihdr['count'] = count
		prihdr['batchN'] = Nb
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		out_filename = outfile_base+str(Nb).zfill(3)+".fits"
		print( out_filename )
		os.system("rm "+out_filename)
		thdulist.writeto(out_filename)
		
	def cornerLCpositionCatalog(self, ii, DMIN=0., DMAX=1000., vmin=190, vmax=100000, NperBatch = 10000000):
		"""
		Extracts the positions and velocity out of a snapshot of the Multidark simulation.   		
		:param ii: index of the snapshot in the list self.snl
		:param DMIN:maximum distance for the pointto be included.
		:param DMAX:maximum distance for the pointto be included.
		:param vmin: name of the quantity of interest, mass, velocity.
		:param vmax: of the quantity of interest in the snapshots.
		:param NperBatch: number of line per fits file, default: 1000000
		 """		
		fl = fileinput.input(self.snl[ii])
		#print self.snl[ii]
		nameSnapshot = self.snl[ii].split('/')[-1][:-5]
		Nb = 0
		count = 0
		output = n.empty((NperBatch,11))
		for line in fl:
			if line[0] == "#" :
				continue

			line = line.split()
			newline =n.array([ float(line[self.columnDict['x']]), float(line[self.columnDict['y']]), float(line[self.columnDict['z']]), float(line[self.columnDict['vx']]), float(line[self.columnDict['vy']]), float(line[self.columnDict['vz']]), float(line[self.columnDict['vmax']]), float(line[self.columnDict['Vpeak']]), n.log10(float(line[self.columnDict['mvir']])), float(line[self.columnDict['rvir']]), float(line[self.columnDict['pid']]) ])
			distance = (newline[0]**2+newline[1]**2+newline[2]**2)**0.5
			if float(line[self.columnDict['vmax']])>vmin and float(line[self.columnDict['vmax']])<vmax and distance<DMAX and distance>=DMIN :
				output[count] = newline
				count+=1
				
			if count == NperBatch  :
				#print "count",count
				#print output
				#print output.shape
				#print output.T[0].shape
				#define the columns
				col0 = fits.Column(name='x',format='D', array=output.T[0] )
				col1 = fits.Column(name='y',format='D', array= output.T[1] )
				col2 = fits.Column(name='z',format='D', array= output.T[2] )
				col3 = fits.Column(name='vx',format='D', array=output.T[3] )
				col4 = fits.Column(name='vy',format='D', array= output.T[4] )
				col5 = fits.Column(name='vz',format='D', array= output.T[5] )
				col6 = fits.Column(name='vmax',format='D', array= output.T[6] )
				col7 = fits.Column(name='vpeak',format='D', array= output.T[7] )
				col8 = fits.Column(name='mvir',format='D', array= output.T[8] )
				col9 = fits.Column(name='rvir',format='D', array= output.T[9] )
				col10 = fits.Column(name='pid',format='D', array= output.T[10] )
				#define the table hdu 
				hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				#define the header
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = nameSnapshot
				prihdr['count'] = count
				prihdr['batchN'] = Nb
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				os.system("rm "+self.snl[ii][:-5]+"_cornerLC_Nb_"+str(Nb)+".fits")
				thdulist.writeto(self.snl[ii][:-5]+"_cornerLC_Nb_"+str(Nb)+".fits")
				Nb+=1
				count=0
				output = n.empty((NperBatch,11))

		# and for the last batch :
		col0 = fits.Column(name='x',format='D', array=output.T[0] )
		col1 = fits.Column(name='y',format='D', array= output.T[1] )
		col2 = fits.Column(name='z',format='D', array= output.T[2] )
		col3 = fits.Column(name='vx',format='D', array=output.T[3] )
		col4 = fits.Column(name='vy',format='D', array= output.T[4] )
		col5 = fits.Column(name='vz',format='D', array= output.T[5] )
		col6 = fits.Column(name='vmax',format='D', array= output.T[6] )
		col7 = fits.Column(name='vpeak',format='D', array= output.T[7] )
		col8 = fits.Column(name='mvir',format='D', array= output.T[8] )
		col9 = fits.Column(name='rvir',format='D', array= output.T[9] )
		col10 = fits.Column(name='pid',format='D', array= output.T[10] )
		#define the table hdu 
		hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		#define the header
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = nameSnapshot
		prihdr['batchN'] = Nb
		prihdr['count'] = count
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		os.system("rm "+self.snl[ii][:-5]+"_cornerLC_Nb_"+str(Nb)+".fits")
		thdulist.writeto(self.snl[ii][:-5]+"_cornerLC_Nb_"+str(Nb)+".fits")
	
	def writeCLUSTERcatalog(self, ii, mmin=10**8, NperBatch = 2000000, file_identifier = "_cluster_Nb_"):
		"""
		Extracts the positions and mass out of a snapshot of the Multidark simulation.        
		:param ii: index of the snapshot in the list self.snl
		:param vmin: name of the quantity of interest, mass, velocity.
		:param vmax: of the quantity of interest in the snapshots.
		:param NperBatch: number of line per fits file, default: 1000000
		 """		
		fl = fileinput.input(self.snl[ii])
		nameSnapshot = os.path.basename(self.snl[ii])[:-5]
		Nb = 0
		count = 0
		output = n.zeros((NperBatch,29))
		for line in fl:
			#print line
			if line[0] == "#" :
				continue

			line = line.split()
			#print line
			#print len(line)
			newline =n.array([int(line[self.columnDict['id']]), float(line[self.columnDict['vmax']]), float(line[self.columnDict['vrms']]), float(line[self.columnDict['rvir']]), float(line[self.columnDict['rs']]), float(line[self.columnDict['pid']]), float(line[self.columnDict['x']]), float(line[self.columnDict['y']]), float(line[self.columnDict['z']]), float(line[self.columnDict['vx']]), float(line[self.columnDict['vy']]), float(line[self.columnDict['vz']]), float(line[self.columnDict['Jx']]), float(line[self.columnDict['Jy']]), float(line[self.columnDict['Jz']]), float(line[self.columnDict['Spin']]), float(line[self.columnDict['Rs_Klypin']]), float(line[self.columnDict['Xoff']]), float(line[self.columnDict['Voff']]), float(line[self.columnDict['Spin_Bullock']]), float(line[self.columnDict['Ax']]), float(line[self.columnDict['Ay']]), float(line[self.columnDict['Az']]), float(line[self.columnDict['b_to_a']]), float(line[self.columnDict['c_to_a']]), n.log10(float(line[self.columnDict['mvir']])), n.log10(float(line[self.columnDict['M200c']])), n.log10(float(line[self.columnDict['M500c']])), n.log10(float(line[self.columnDict['M2500c']])) ])
			#print newline
			#print newline.shape
			if float(line[self.columnDict['M500c']])>mmin :
				output[count] = newline
				count+=1
				
			if count == NperBatch  :
				#print "count",count
				#print output
				#print output.shape
				#print output.T[0].shape
				#define the columns
				hdu_cols  = fits.ColDefs([
				fits.Column(name='id',format='I',            array= output.T[0] )
				,fits.Column(name='vmax',format='D',         array= output.T[1] ) 
				,fits.Column(name='vrms',format='D',         array= output.T[2] ) 
				,fits.Column(name='rvir',format='D',         array= output.T[3] ) 
				,fits.Column(name='rs',format='D',           array= output.T[4] ) 
				,fits.Column(name='pid',format='D',          array= output.T[5] ) 
				,fits.Column(name='x',format='D',            array= output.T[6] ) 
				,fits.Column(name='y',format='D',            array= output.T[7] ) 
				,fits.Column(name='z',format='D',            array= output.T[8] ) 
				,fits.Column(name='vx',format='D',           array= output.T[9] ) 
				,fits.Column(name='vy',format='D',           array= output.T[10] ) 
				,fits.Column(name='vz',format='D',           array= output.T[11] ) 
				,fits.Column(name='Jx',format='D',           array= output.T[12] ) 
				,fits.Column(name='Jy',format='D',           array= output.T[13] ) 
				,fits.Column(name='Jz',format='D',           array= output.T[14] ) 
				,fits.Column(name='Spin',format='D',         array= output.T[15] ) 
				,fits.Column(name='Rs_Klypin',format='D',    array= output.T[16] ) 
				,fits.Column(name='Xoff',format='D',         array= output.T[17] ) 
				,fits.Column(name='Voff',format='D',         array= output.T[18] ) 
				,fits.Column(name='Spin_Bullock',format='D', array= output.T[19] ) 
				,fits.Column(name='Ax',format='D',           array= output.T[20] ) 
				,fits.Column(name='Ay',format='D',           array= output.T[21] ) 
				,fits.Column(name='Az',format='D',           array= output.T[22] ) 
				,fits.Column(name='b_to_a',format='D',       array= output.T[23] ) 
				,fits.Column(name='c_to_a',format='D',       array= output.T[24] ) 
				,fits.Column(name='mvir',format='D',         array= output.T[25] ) 
				,fits.Column(name='M200c',format='D',        array= output.T[26] ) 
				,fits.Column(name='M500c',format='D',        array= output.T[27] ) 
				,fits.Column(name='M2500c',format='D',       array= output.T[28] )
				])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				#define the header
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = nameSnapshot
				prihdr['count'] = count
				prihdr['batchN'] = Nb
				prihdr['author'] = 'JC'
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				os.system("rm "+self.snl[ii][:-5]+file_identifier+str(Nb)+".fits")
				thdulist.writeto(self.snl[ii][:-5]+file_identifier+str(Nb)+".fits")
				Nb+=1
				count=0
				#resest the output matrix
				output = n.zeros((NperBatch,29))
		
		
		# and for the last batch :		
		hdu_cols  = fits.ColDefs([
		fits.Column(name='id',format='I',            array= output.T[0] )
		,fits.Column(name='vmax',format='D',         array= output.T[1] ) 
		,fits.Column(name='vrms',format='D',         array= output.T[2] ) 
		,fits.Column(name='rvir',format='D',         array= output.T[3] ) 
		,fits.Column(name='rs',format='D',           array= output.T[4] ) 
		,fits.Column(name='pid',format='D',          array= output.T[5] ) 
		,fits.Column(name='x',format='D',            array= output.T[6] ) 
		,fits.Column(name='y',format='D',            array= output.T[7] ) 
		,fits.Column(name='z',format='D',            array= output.T[8] ) 
		,fits.Column(name='vx',format='D',           array= output.T[9] ) 
		,fits.Column(name='vy',format='D',           array= output.T[10] ) 
		,fits.Column(name='vz',format='D',           array= output.T[11] ) 
		,fits.Column(name='Jx',format='D',           array= output.T[12] ) 
		,fits.Column(name='Jy',format='D',           array= output.T[13] ) 
		,fits.Column(name='Jz',format='D',           array= output.T[14] ) 
		,fits.Column(name='Spin',format='D',         array= output.T[15] ) 
		,fits.Column(name='Rs_Klypin',format='D',    array= output.T[16] ) 
		,fits.Column(name='Xoff',format='D',         array= output.T[17] ) 
		,fits.Column(name='Voff',format='D',         array= output.T[18] ) 
		,fits.Column(name='Spin_Bullock',format='D', array= output.T[19] ) 
		,fits.Column(name='Ax',format='D',           array= output.T[20] ) 
		,fits.Column(name='Ay',format='D',           array= output.T[21] ) 
		,fits.Column(name='Az',format='D',           array= output.T[22] ) 
		,fits.Column(name='b_to_a',format='D',       array= output.T[23] ) 
		,fits.Column(name='c_to_a',format='D',       array= output.T[24] ) 
		,fits.Column(name='mvir',format='D',         array= output.T[25] ) 
		,fits.Column(name='M200c',format='D',        array= output.T[26] ) 
		,fits.Column(name='M500c',format='D',        array= output.T[27] ) 
		,fits.Column(name='M2500c',format='D',       array= output.T[28] )
		])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		#define the header
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = nameSnapshot
		prihdr['count'] = count
		prihdr['batchN'] = Nb
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		os.system("rm "+self.snl[ii][:-5]+file_identifier+str(Nb)+".fits")
		thdulist.writeto(self.snl[ii][:-5]+file_identifier+str(Nb)+".fits")
	
	def writeSAMcatalog(self, ii, mmin=10**8, NperBatch = 2000000, file_identifier = "_SAM_Nb_"):
		"""
		Extracts the positions and mass out of a snapshot of the Multidark simulation.        
		:param ii: index of the snapshot in the list self.snl
		:param vmin: name of the quantity of interest, mass, velocity.
		:param vmax: of the quantity of interest in the snapshots.
		:param NperBatch: number of line per fits file, default: 1000000
		 """		
		fl = fileinput.input(self.snl[ii])
		nameSnapshot = os.path.basename(self.snl[ii])[:-5]
		Nb = 0
		count = 0
		output = n.zeros((NperBatch,29))
		for line in fl:
			#print line
			if line[0] == "#" :
				continue

			line = line.split()
			#print line
			#print len(line)
			newline =n.array([int(line[self.columnDict['id']]), float(line[self.columnDict['vmax']]), float(line[self.columnDict['vrms']]), float(line[self.columnDict['rvir']]), float(line[self.columnDict['rs']]), float(line[self.columnDict['pid']]), float(line[self.columnDict['x']]), float(line[self.columnDict['y']]), float(line[self.columnDict['z']]), float(line[self.columnDict['vx']]), float(line[self.columnDict['vy']]), float(line[self.columnDict['vz']]), float(line[self.columnDict['Jx']]), float(line[self.columnDict['Jy']]), float(line[self.columnDict['Jz']]), float(line[self.columnDict['Spin']]), float(line[self.columnDict['Rs_Klypin']]), float(line[self.columnDict['Xoff']]), float(line[self.columnDict['Voff']]), float(line[self.columnDict['Spin_Bullock']]), float(line[self.columnDict['Ax']]), float(line[self.columnDict['Ay']]), float(line[self.columnDict['Az']]), float(line[self.columnDict['b_to_a']]), float(line[self.columnDict['c_to_a']]), n.log10(float(line[self.columnDict['mvir']])), n.log10(float(line[self.columnDict['M200c']])), n.log10(float(line[self.columnDict['M500c']])), n.log10(float(line[self.columnDict['M2500c']])) ])
			#print newline
			#print newline.shape
			if float(line[self.columnDict['mvir']])>mmin :
				output[count] = newline
				count+=1
				
			if count == NperBatch  :
				out_filename = os.path.join(os.environ["MD10"],"work_agn",os.path.basename(self.snl[ii])[:-5]+file_identifier+str(Nb)+".fits")
				print( out_filename )
				#print "count",count
				#print output
				#print output.shape
				#print output.T[0].shape
				#define the columns
				hdu_cols  = fits.ColDefs([
				fits.Column(name='id',format='I',            array= output.T[0] )
				,fits.Column(name='vmax',format='D',         array= output.T[1] ) 
				,fits.Column(name='vrms',format='D',         array= output.T[2] ) 
				,fits.Column(name='rvir',format='D',         array= output.T[3] ) 
				,fits.Column(name='rs',format='D',           array= output.T[4] ) 
				,fits.Column(name='pid',format='D',          array= output.T[5] ) 
				,fits.Column(name='x',format='D',            array= output.T[6] ) 
				,fits.Column(name='y',format='D',            array= output.T[7] ) 
				,fits.Column(name='z',format='D',            array= output.T[8] ) 
				,fits.Column(name='vx',format='D',           array= output.T[9] ) 
				,fits.Column(name='vy',format='D',           array= output.T[10] ) 
				,fits.Column(name='vz',format='D',           array= output.T[11] ) 
				,fits.Column(name='Jx',format='D',           array= output.T[12] ) 
				,fits.Column(name='Jy',format='D',           array= output.T[13] ) 
				,fits.Column(name='Jz',format='D',           array= output.T[14] ) 
				,fits.Column(name='Spin',format='D',         array= output.T[15] ) 
				,fits.Column(name='Rs_Klypin',format='D',    array= output.T[16] ) 
				,fits.Column(name='Xoff',format='D',         array= output.T[17] ) 
				,fits.Column(name='Voff',format='D',         array= output.T[18] ) 
				,fits.Column(name='Spin_Bullock',format='D', array= output.T[19] ) 
				,fits.Column(name='Ax',format='D',           array= output.T[20] ) 
				,fits.Column(name='Ay',format='D',           array= output.T[21] ) 
				,fits.Column(name='Az',format='D',           array= output.T[22] ) 
				,fits.Column(name='b_to_a',format='D',       array= output.T[23] ) 
				,fits.Column(name='c_to_a',format='D',       array= output.T[24] ) 
				,fits.Column(name='mvir',format='D',         array= output.T[25] ) 
				,fits.Column(name='M200c',format='D',        array= output.T[26] ) 
				,fits.Column(name='M500c',format='D',        array= output.T[27] ) 
				,fits.Column(name='M2500c',format='D',       array= output.T[28] )
				])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				#define the header
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = nameSnapshot
				prihdr['count'] = count
				prihdr['batchN'] = Nb
				prihdr['author'] = 'JC'
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				out_filename = os.path.join(os.environ["MD10"],"work_agn",os.path.basename(self.snl[ii])[:-5]+file_identifier+str(Nb)+".fits")
				print( out_filename )
				os.system("rm "+out_filename)
				thdulist.writeto(out_filename)
				Nb+=1
				count=0
				#resest the output matrix
				output = n.zeros((NperBatch,29))
		
		
		# and for the last batch :		
		hdu_cols  = fits.ColDefs([
		fits.Column(name='id',format='I',            array= output.T[0] )
		,fits.Column(name='vmax',format='D',         array= output.T[1] ) 
		,fits.Column(name='vrms',format='D',         array= output.T[2] ) 
		,fits.Column(name='rvir',format='D',         array= output.T[3] ) 
		,fits.Column(name='rs',format='D',           array= output.T[4] ) 
		,fits.Column(name='pid',format='D',          array= output.T[5] ) 
		,fits.Column(name='x',format='D',            array= output.T[6] ) 
		,fits.Column(name='y',format='D',            array= output.T[7] ) 
		,fits.Column(name='z',format='D',            array= output.T[8] ) 
		,fits.Column(name='vx',format='D',           array= output.T[9] ) 
		,fits.Column(name='vy',format='D',           array= output.T[10] ) 
		,fits.Column(name='vz',format='D',           array= output.T[11] ) 
		,fits.Column(name='Jx',format='D',           array= output.T[12] ) 
		,fits.Column(name='Jy',format='D',           array= output.T[13] ) 
		,fits.Column(name='Jz',format='D',           array= output.T[14] ) 
		,fits.Column(name='Spin',format='D',         array= output.T[15] ) 
		,fits.Column(name='Rs_Klypin',format='D',    array= output.T[16] ) 
		,fits.Column(name='Xoff',format='D',         array= output.T[17] ) 
		,fits.Column(name='Voff',format='D',         array= output.T[18] ) 
		,fits.Column(name='Spin_Bullock',format='D', array= output.T[19] ) 
		,fits.Column(name='Ax',format='D',           array= output.T[20] ) 
		,fits.Column(name='Ay',format='D',           array= output.T[21] ) 
		,fits.Column(name='Az',format='D',           array= output.T[22] ) 
		,fits.Column(name='b_to_a',format='D',       array= output.T[23] ) 
		,fits.Column(name='c_to_a',format='D',       array= output.T[24] ) 
		,fits.Column(name='mvir',format='D',         array= output.T[25] ) 
		,fits.Column(name='M200c',format='D',        array= output.T[26] ) 
		,fits.Column(name='M500c',format='D',        array= output.T[27] ) 
		,fits.Column(name='M2500c',format='D',       array= output.T[28] )
		])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		#define the header
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = nameSnapshot
		prihdr['count'] = count
		prihdr['batchN'] = Nb
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		out_filename = os.path.join(os.environ["MD10"],"work_agn",os.path.basename(self.snl[ii])[:-5]+file_identifier+str(Nb)+".fits")
		print( out_filename )
		os.system("rm "+out_filename)
		thdulist.writeto(out_filename)
		
	def writePositionCatalogPM(self, ii, vmin=30., mmin=10**8, NperBatch = 20000000):
		"""
		Extracts the positions and velocity out of a snapshot of the Multidark simulation.        
		:param ii: index of the snapshot in the list self.snl
		:param vmin: name of the quantity of interest, mass, velocity.
		:param vmax: of the quantity of interest in the snapshots.
		:param NperBatch: number of line per fits file, default: 1000000
		 """		
		fl = fileinput.input(self.snl[ii])
		nameSnapshot = self.snl[ii].split('/')[-1][:-5]
		Nb = 0
		count = 0
		output = n.zeros((NperBatch,7))
		for line in fl:
			if line[0] == "#" :
				continue

			line = line.split()
			newline =n.array([int(line[self.columnDict['id']]), float(line[self.columnDict['pid']]), float(line[self.columnDict['x']]), float(line[self.columnDict['y']]), float(line[self.columnDict['z']]), float(line[self.columnDict['vmax']]), n.log10(float(line[self.columnDict['mvir']])) ])
			if   float(line[self.columnDict['vmax']])>vmin and float(line[self.columnDict['mvir']])>mmin :
				output[count] = newline
				count+=1
				
			if count == NperBatch  :
				#print "count",count
				#print output
				#print output.shape
				#print output.T[0].shape
				#define the columns
				col0 = fits.Column(name='id',format='D', array= output.T[0] )
				col1 = fits.Column(name='pid',format='D', array= output.T[1] )
				col2 = fits.Column(name='x',format='D', array=output.T[2] )
				col3 = fits.Column(name='y',format='D', array= output.T[3] )
				col4 = fits.Column(name='z',format='D', array= output.T[4] )
				col5 = fits.Column(name='vmax',format='D', array= output.T[5] )
				col6 = fits.Column(name='mvir',format='D', array=output.T[6] )
				#define the table hdu 
				hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				#define the header
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = nameSnapshot
				prihdr['count'] = count
				prihdr['batchN'] = Nb
				prihdr['author'] = 'JC'
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				os.system("rm "+self.snl[ii][:-5]+"_PM_Nb_"+str(Nb)+".fits")
				thdulist.writeto(self.snl[ii][:-5]+"_PM_Nb_"+str(Nb)+".fits")
				Nb+=1
				count=0
				#resest the output matrix
				output = n.zeros((NperBatch,7))
		
		
		# and for the last batch :		
		col0 = fits.Column(name='id',format='D', array= output.T[0][:count] )
		col1 = fits.Column(name='pid',format='D', array= output.T[1][:count] )
		col2 = fits.Column(name='x',format='D', array=output.T[2][:count] )
		col3 = fits.Column(name='y',format='D', array= output.T[3][:count] )
		col4 = fits.Column(name='z',format='D', array= output.T[4][:count] )
		col5 = fits.Column(name='vmax',format='D', array= output.T[5][:count] )
		col6 = fits.Column(name='mvir',format='D', array=output.T[6][:count] )
		#define the table hdu 
		hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5, col6])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		#define the header
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = nameSnapshot
		prihdr['count'] = count
		prihdr['batchN'] = Nb
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		os.system("rm "+self.snl[ii][:-5]+"_PM_Nb_"+str(Nb)+".fits")
		thdulist.writeto(self.snl[ii][:-5]+"_PM_Nb_"+str(Nb)+".fits")
	
	def writePositionCatalogVmaxM200c(self, ii, vmin=190, vmax=10000, NperBatch = 10000000):
		"""
		Extracts the positions and velocity out of a snapshot of the Multidark simulation.        
		:param ii: index of the snapshot in the list self.snl
		:param vmin: name of the quantity of interest, mass, velocity.
		:param vmax: of the quantity of interest in the snapshots.
		:param NperBatch: number of line per fits file, default: 1000000
		 """		
		fl = fileinput.input(self.snl[ii])
		nameSnapshot = self.snl[ii].split('/')[-1][:-5]
		Nb = 0
		count = 0
		output = n.empty((NperBatch,6))
		for line in fl:
			if line[0] == "#" :
				continue

			line = line.split()
			newline =n.array([ float(line[self.columnDict['x']]), float(line[self.columnDict['y']]), float(line[self.columnDict['z']]), float(line[self.columnDict['vmax']]), n.log10(float(line[self.columnDict['M200c']])), float(line[self.columnDict['pid']]) ])
			if newline[3]>vmin and newline[3]<vmax :
				output[count] = newline
				count+=1
				
			if count == NperBatch  :
				#print "count",count
				#print output
				#print output.shape
				#print output.T[0].shape
				#define the columns
				col0 = fits.Column(name='x',format='D', array=output.T[0] )
				col1 = fits.Column(name='y',format='D', array= output.T[1] )
				col2 = fits.Column(name='z',format='D', array= output.T[2] )
				col3 = fits.Column(name='vmax',format='D', array= output.T[3] )
				col4 = fits.Column(name='M200c',format='D', array= output.T[4] )
				col5 = fits.Column(name='pid',format='D', array= output.T[5] )
				#define the table hdu 
				hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				#define the header
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = nameSnapshot
				prihdr['count'] = count
				prihdr['batchN'] = Nb
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				os.system("rm "+self.snl[ii][:-5]+"_VmaxM200c_Nb_"+str(Nb)+".fits")
				thdulist.writeto(self.snl[ii][:-5]+"_VmaxM200c_Nb_"+str(Nb)+".fits")
				Nb+=1
				count=0
				output = n.empty((NperBatch,6))

		# and for the last batch :
		col0 = fits.Column(name='x',format='D', array=output.T[0] )
		col1 = fits.Column(name='y',format='D', array= output.T[1] )
		col2 = fits.Column(name='z',format='D', array= output.T[2] )
		col3 = fits.Column(name='vmax',format='D', array= output.T[3] )
		col4 = fits.Column(name='M200c',format='D', array= output.T[4] )
		col5 = fits.Column(name='pid',format='D', array= output.T[5] )
		#define the table hdu 
		hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4, col5])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		#define the header
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = nameSnapshot
		prihdr['count'] = count
		prihdr['batchN'] = Nb
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		os.system("rm "+self.snl[ii][:-5]+"_VmaxM200c_Nb_"+str(Nb)+".fits")
		thdulist.writeto(self.snl[ii][:-5]+"_VmaxM200c_Nb_"+str(Nb)+".fits")
	
	def writePositionCatalogVmax(self, ii, vmin=190, vmax=10000, NperBatch = 10000000):
		"""
		Extracts the positions and velocity out of a snapshot of the Multidark simulation.        
		:param ii: index of the snapshot in the list self.snl
		:param vmin: name of the quantity of interest, mass, velocity.
		:param vmax: of the quantity of interest in the snapshots.
		:param NperBatch: number of line per fits file, default: 1000000
		 """		
		fl = fileinput.input(self.snl[ii])
		nameSnapshot = self.snl[ii].split('/')[-1][:-5]
		Nb = 0
		count = 0
		output = n.empty((NperBatch,5))
		for line in fl:
			if line[0] == "#" :
				continue

			line = line.split()
			newline =n.array([ float(line[self.columnDict['x']]), float(line[self.columnDict['y']]), float(line[self.columnDict['z']]), float(line[self.columnDict['vmax']]), float(line[self.columnDict['pid']]) ])
			if newline[3]>vmin and newline[3]<vmax :
				output[count] = newline
				count+=1
				
			if count == NperBatch  :
				#print "count",count
				#print output
				#print output.shape
				#print output.T[0].shape
				#define the columns
				col0 = fits.Column(name='x',format='D', array=output.T[0] )
				col1 = fits.Column(name='y',format='D', array= output.T[1] )
				col2 = fits.Column(name='z',format='D', array= output.T[2] )
				col3 = fits.Column(name='vmax',format='D', array= output.T[3] )
				col4 = fits.Column(name='pid',format='D', array= output.T[4] )
				#define the table hdu 
				hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4])
				tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
				#define the header
				prihdr = fits.Header()
				prihdr['HIERARCH nameSnapshot'] = nameSnapshot
				prihdr['count'] = count
				prihdr['batchN'] = Nb
				prihdu = fits.PrimaryHDU(header=prihdr)
				#writes the file
				thdulist = fits.HDUList([prihdu, tb_hdu])
				os.system("rm "+self.snl[ii][:-5]+"_Nb_"+str(Nb)+".fits")
				thdulist.writeto(self.snl[ii][:-5]+"_Nb_"+str(Nb)+".fits")
				Nb+=1
				count=0
				output = n.empty((NperBatch,5))

		# and for the last batch :
		col0 = fits.Column(name='x',format='D', array= output.T[0][:count])
		col1 = fits.Column(name='y',format='D', array= output.T[1][:count])
		col2 = fits.Column(name='z',format='D', array= output.T[2][:count])
		col3 = fits.Column(name='vmax',format='D', array= output.T[3][:count])
		col4 = fits.Column(name='pid',format='D', array= output.T[4][:count])
		#define the table hdu 
		hdu_cols  = fits.ColDefs([col0, col1, col2, col3, col4])
		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		#define the header
		prihdr = fits.Header()
		prihdr['HIERARCH nameSnapshot'] = nameSnapshot
		prihdr['batchN'] = Nb
		prihdr['count'] = count
		prihdu = fits.PrimaryHDU(header=prihdr)
		#writes the file
		thdulist = fits.HDUList([prihdu, tb_hdu])
		os.system("rm "+self.snl[ii][:-5]+"_Nb_"+str(Nb)+".fits")
		thdulist.writeto(self.snl[ii][:-5]+"_Nb_"+str(Nb)+".fits")
	
	def compute2PCF(self, catalogList, vmin=65, rmax=200, dlogBin=0.05, Nmax=4000000., dr = 1., name = ""):
		"""
		Extracts the 2PCF out of a catalog of halos        
		:param catalog: where the catalog is
		:param vmin: minimum circular velocity.
		:param dlogBin: bin width.
		:param rmax: maximum distance
		"""
		hdus = []
		for ii in n.arange(len(catalogList)):
			hdus.append( fits.open(catalogList[ii])[1].data )

		vbins = 10**n.arange(n.log10(vmin),4. ,dlogBin)
		for jj in range(len(vbins)-1):
			outfile = catalogList[0][:-10] + "_vmax_" +str(n.round(vbins[jj],2))+ "_" +str(n.round(vbins[jj+1],2)) + "_" + name + "_xiR.pkl"
			t0 = time.time()
			sel = n.array([ (hdu['vmax']>vbins[jj])&(hdu['vmax']<vbins[jj+1]) for hdu in hdus])
			xR = n.hstack(( n.array([ hdus[ii]['x'][sel[ii]] for ii in range(len(hdus)) ]) ))
			yR = n.hstack(( n.array([ hdus[ii]['y'][sel[ii]] for ii in range(len(hdus)) ]) ))
			zR = n.hstack(( n.array([ hdus[ii]['z'][sel[ii]] for ii in range(len(hdus)) ]) ))
			Ntotal = len(xR)
			if len(xR)>5000 and len(xR)<=Nmax:
				#print vbins[jj], vbins[jj+1]
				insideSel=(xR>rmax)&(xR<self.Lbox.value-rmax)&(yR>rmax)&(yR<self.Lbox.value-rmax)&(zR>rmax)&(zR<self.Lbox.value-rmax)
				volume=(self.Lbox.value)**3
				# defines the trees
				#print "creates trees"
				treeRandoms=t.cKDTree(n.transpose([xR,yR,zR]),1000.0)
				treeData=t.cKDTree(n.transpose([xR[insideSel],yR[insideSel],zR[insideSel]]),1000.0)
				nD=len(treeData.data)
				nR=len(treeRandoms.data)
				#print nD, nR
				bin_xi3D=n.arange(0, rmax, dr)
				# now does the pair counts :
				pairs=treeData.count_neighbors(treeRandoms, bin_xi3D)
				t3 = time.time()
				DR=pairs[1:]-pairs[:-1]
				dV= (bin_xi3D[1:]**3 - bin_xi3D[:-1]**3 )*4*n.pi/3.
				pairCount=nD*nR#-nD*(nD-1)/2.
				xis = DR*volume/(dV * pairCount) -1.
				f=open(outfile,'w')
				cPickle.dump([bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbins[jj], vbins[jj+1]],f)
				f.close()
				t4 = time.time()
				#print "total time in s, min",t4 - t0, (t4 - t0)/60.
				#return DR, volume, dV, pairCount, pairs, nD, nR

			if  len(xR)>Nmax:
				#print vbins[jj], vbins[jj+1], "downsampling ..."
				downSamp = (n.random.random(len(xR))<Nmax / float(len(xR)) )
				xR = xR[downSamp]
				yR = yR[downSamp]
				zR = zR[downSamp]
				
				insideSel=(xR>rmax)&(xR<self.Lbox.value-rmax)&(yR>rmax)&(yR<self.Lbox.value-rmax)&(zR>rmax)&(zR<self.Lbox.value-rmax)
				volume=(self.Lbox.value-rmax*2)**3
				# defines the trees
				#print "creates trees"
				treeRandoms=t.cKDTree(n.transpose([xR,yR,zR]),1000.0)
				treeData=t.cKDTree(n.transpose([xR[insideSel],yR[insideSel],zR[insideSel]]),1000.0)
				nD=len(treeData.data)
				nR=len(treeRandoms.data)
				#print nD, nR
				bin_xi3D=n.arange(0, rmax, dr)
				# now does the pair counts :
				pairs=treeData.count_neighbors(treeRandoms, bin_xi3D)
				t3 = time.time()
				DR=pairs[1:]-pairs[:-1]
				dV=4*n.pi*(bin_xi3D[1:]**3 - bin_xi3D[:-1]**3 )/3.
				pairCount=nD*nR#-nD*(nD-1)/2.
				xis = DR*volume/(dV * pairCount) -1.
				f=open(outfile,'w')
				cPickle.dump([bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbins[jj], vbins[jj+1]],f)
				f.close()
				t4 = time.time()
				#print "total time in s, min",t4 - t0, (t4 - t0)/60.
				#return DR, volume, dV, pairCount, pairs, nD, nR

	def compute2PCF_MASS(self, catalogList, rmax=200, dr = 0.1, vmin=9, dlogBin=0.05, Nmax=2000000.,  name = ""):
		"""
		Extracts the 2PCF out of a catalog of halos        
		:param catalog: where the catalog is
		:param vmin: minimum circular velocity.
		:param dlogBin: bin width.
		:param rmax: maximum distance
		"""
		hdus = []
		for ii in n.arange(len(catalogList)):
			hdus.append( fits.open(catalogList[ii])[1].data )

		vbins = n.arange(vmin, 16. ,dlogBin)
		#n.arange(8,16,0.05)
		for jj in range(len(vbins)-1):
			outfile = catalogList[0][:-13] + "_mvir_" +str(n.round(vbins[jj],2))+ "_" +str(n.round(vbins[jj+1],2)) + "_" + name + "_xiR.pkl"
			print outfile
			if os.path.isfile(outfile):
				print "pass"
				pass
			else :
				t0 = time.time()
				sel = n.array([ (hdu['mvir']>vbins[jj])&(hdu['mvir']<vbins[jj+1]) for hdu in hdus])
				xR = n.hstack(( n.array([ hdus[ii]['x'][sel[ii]] for ii in range(len(hdus)) ]) ))
				yR = n.hstack(( n.array([ hdus[ii]['y'][sel[ii]] for ii in range(len(hdus)) ]) ))
				zR = n.hstack(( n.array([ hdus[ii]['z'][sel[ii]] for ii in range(len(hdus)) ]) ))
				Ntotal = len(xR)
				if len(xR)>10000 and len(xR)<=Nmax:
					#print vbins[jj], vbins[jj+1]
					insideSel=(xR>rmax)&(xR<self.Lbox.value-rmax)&(yR>rmax)&(yR<self.Lbox.value-rmax)&(zR>rmax)&(zR<self.Lbox.value-rmax)
					volume=(self.Lbox.value)**3
					# defines the trees
					#print "creates trees"
					treeRandoms=t.cKDTree(n.transpose([xR,yR,zR]),1000.0)
					treeData=t.cKDTree(n.transpose([xR[insideSel],yR[insideSel],zR[insideSel]]),1000.0)
					nD=len(treeData.data)
					nR=len(treeRandoms.data)
					#print nD, nR
					bin_xi3D=n.arange(0, rmax, dr)
					# now does the pair counts :
					pairs=treeData.count_neighbors(treeRandoms, bin_xi3D)
					t3 = time.time()
					DR=pairs[1:]-pairs[:-1]
					dV= (bin_xi3D[1:]**3 - bin_xi3D[:-1]**3 )*4*n.pi/3.
					pairCount=nD*nR#-nD*(nD-1)/2.
					xis = DR*volume/(dV * pairCount) -1.
					f=open(outfile,'w')
					cPickle.dump([bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbins[jj], vbins[jj+1]],f)
					f.close()
					t4 = time.time()
					#print "total time in s, min",t4 - t0, (t4 - t0)/60.
					#return DR, volume, dV, pairCount, pairs, nD, nR

				if  len(xR)>Nmax:
					#print vbins[jj], vbins[jj+1], "downsampling ..."
					downSamp = (n.random.random(len(xR))<Nmax / float(len(xR)) )
					xR = xR[downSamp]
					yR = yR[downSamp]
					zR = zR[downSamp]
					
					insideSel=(xR>rmax)&(xR<self.Lbox.value-rmax)&(yR>rmax)&(yR<self.Lbox.value-rmax)&(zR>rmax)&(zR<self.Lbox.value-rmax)
					volume=(self.Lbox.value-rmax*2)**3
					# defines the trees
					#print "creates trees"
					treeRandoms=t.cKDTree(n.transpose([xR,yR,zR]),1000.0)
					treeData=t.cKDTree(n.transpose([xR[insideSel],yR[insideSel],zR[insideSel]]),1000.0)
					nD=len(treeData.data)
					nR=len(treeRandoms.data)
					#print nD, nR
					bin_xi3D=n.arange(0, rmax, dr)
					# now does the pair counts :
					pairs=treeData.count_neighbors(treeRandoms, bin_xi3D)
					t3 = time.time()
					DR=pairs[1:]-pairs[:-1]
					dV=4*n.pi*(bin_xi3D[1:]**3 - bin_xi3D[:-1]**3 )/3.
					pairCount=nD*nR#-nD*(nD-1)/2.
					xis = DR*volume/(dV * pairCount) -1.
					f=open(outfile,'w')
					cPickle.dump([bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbins[jj], vbins[jj+1]],f)
					f.close()
					t4 = time.time()
					#print "total time in s, min",t4 - t0, (t4 - t0)/60.
					#return DR, volume, dV, pairCount, pairs, nD, nR

	def compute2PCF_LX(self, catalogList, rmax=200, dr = 0.1, vmin=41., dlogBin=0.25, Nmax=2000000.,  name = ""):
		"""
		Extracts the 2PCF out of a catalog of halos        
		:param catalog: where the catalog is
		:param vmin: minimum circular velocity.
		:param dlogBin: bin width.
		:param rmax: maximum distance
		"""
		# first create the catalog
		hdus = []
		for ii in n.arange(len(catalogList)):
			hd = fits.open(catalogList[ii])[1].data
			cut_min = (hd['lambda_sar_Bo16']+hd['Mgal_mvir_Mo13'] > vmin)
			hdus.append( hd[cut_min] )
			print catalogList[ii], len(hd[cut_min])

		vbins = n.arange(vmin, 43.51 ,dlogBin)
		#n.arange(8,16,0.05)
		for jj in range(len(vbins)-1):
			outfile = catalogList[0][:-13] + "_LX_" +str(n.round(vbins[jj],2))+ "_" +str(n.round(vbins[jj+1],2)) + "_" + name + "_xiR.pkl"
			print outfile
			if os.path.isfile(outfile):
				print "pass"
				pass
			else :
				t0 = time.time()
				sel = n.array([ (hdu['lambda_sar_Bo16']+hdu['Mgal_mvir_Mo13']>vbins[jj])&(hdu['lambda_sar_Bo16']+hdu['Mgal_mvir_Mo13']<vbins[jj+1]) for hdu in hdus])
				xR = n.hstack(( n.array([ hdus[ii]['x'][sel[ii]] for ii in range(len(hdus)) ]) ))
				yR = n.hstack(( n.array([ hdus[ii]['y'][sel[ii]] for ii in range(len(hdus)) ]) ))
				zR = n.hstack(( n.array([ hdus[ii]['z'][sel[ii]] for ii in range(len(hdus)) ]) ))
				Ntotal = len(xR)
				if len(xR)>10000 and len(xR)<=Nmax:
					#print vbins[jj], vbins[jj+1]
					insideSel=(xR>rmax)&(xR<self.Lbox.value-rmax)&(yR>rmax)&(yR<self.Lbox.value-rmax)&(zR>rmax)&(zR<self.Lbox.value-rmax)
					volume=(self.Lbox.value)**3
					# defines the trees
					#print "creates trees"
					treeRandoms=t.cKDTree(n.transpose([xR,yR,zR]),1000.0)
					treeData=t.cKDTree(n.transpose([xR[insideSel],yR[insideSel],zR[insideSel]]),1000.0)
					nD=len(treeData.data)
					nR=len(treeRandoms.data)
					#print nD, nR
					bin_xi3D=n.arange(0, rmax, dr)
					# now does the pair counts :
					pairs=treeData.count_neighbors(treeRandoms, bin_xi3D)
					t3 = time.time()
					DR=pairs[1:]-pairs[:-1]
					dV= (bin_xi3D[1:]**3 - bin_xi3D[:-1]**3 )*4*n.pi/3.
					pairCount=nD*nR#-nD*(nD-1)/2.
					xis = DR*volume/(dV * pairCount) -1.
					f=open(outfile,'w')
					cPickle.dump([bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbins[jj], vbins[jj+1]],f)
					f.close()
					t4 = time.time()
					#print "total time in s, min",t4 - t0, (t4 - t0)/60.
					#return DR, volume, dV, pairCount, pairs, nD, nR

				if  len(xR)>Nmax:
					#print vbins[jj], vbins[jj+1], "downsampling ..."
					downSamp = (n.random.random(len(xR))<Nmax / float(len(xR)) )
					xR = xR[downSamp]
					yR = yR[downSamp]
					zR = zR[downSamp]
					
					insideSel=(xR>rmax)&(xR<self.Lbox.value-rmax)&(yR>rmax)&(yR<self.Lbox.value-rmax)&(zR>rmax)&(zR<self.Lbox.value-rmax)
					volume=(self.Lbox.value-rmax*2)**3
					# defines the trees
					#print "creates trees"
					treeRandoms=t.cKDTree(n.transpose([xR,yR,zR]),1000.0)
					treeData=t.cKDTree(n.transpose([xR[insideSel],yR[insideSel],zR[insideSel]]),1000.0)
					nD=len(treeData.data)
					nR=len(treeRandoms.data)
					#print nD, nR
					bin_xi3D=n.arange(0, rmax, dr)
					# now does the pair counts :
					pairs=treeData.count_neighbors(treeRandoms, bin_xi3D)
					t3 = time.time()
					DR=pairs[1:]-pairs[:-1]
					dV=4*n.pi*(bin_xi3D[1:]**3 - bin_xi3D[:-1]**3 )/3.
					pairCount=nD*nR#-nD*(nD-1)/2.
					xis = DR*volume/(dV * pairCount) -1.
					f=open(outfile,'w')
					cPickle.dump([bin_xi3D,xis, DR, volume, dV, pairCount, pairs, Ntotal, nD, nR, vbins[jj], vbins[jj+1]],f)
					f.close()
					t4 = time.time()
					#print "total time in s, min",t4 - t0, (t4 - t0)/60.
					#return DR, volume, dV, pairCount, pairs, nD, nR

	def computeSingleDistributionFunctionJKresampling(self, fileList, rootname, name, bins, Ljk = 100., overlap = 1. ) :
		"""
		Extracts the distribution of quantity 'name' out of all snapshots of the Multidark simulation.
		Resamples the box in smaller boxes of length Ljk in Mpc/h
		:param ii: index of the snapshot in the list self.snl
		:param name: name of the quantity of interest, mass, velocity.
		:param index: of the quantity of interest in the snapshots.
		:param bins: binning scheme to compute the historgram.
		:param Ljk: length of the resampled box
		:param overlap: allowed overlap between resampled realizations : 1 = no overlap 2 : 50% overlap ... 
		"""		
		output_dir = join(self.wdir,self.boxDir,"properties",name)
		os.system('mkdir '+ output_dir)
		# define boundaries
		NBoundariesPerSide = int(overlap*self.Lbox.value/Ljk)
		bounds = n.arange(NBoundariesPerSide+1)* Ljk / overlap
		#print "boundaries on each side: ", bounds
		Xi, Yi, Zi = n.meshgrid(bounds[:-1],bounds[:-1],bounds[:-1])
		X = n.ravel(Xi)
		Y = n.ravel(Yi)
		Z = n.ravel(Zi)	
		#print X.min(), X.max(), len(X),len(bounds)
		# loops over the fileList : fits files with the data
		nnC = n.zeros((len(fileList),len(X),len(bins)-1))
		nnS = n.zeros((len(fileList),len(X),len(bins)-1))
		for jj, file in enumerate(fileList):
			#print file
			dd = fits.open(file)[1].data
			cen = (dd['pid']==-1)
			sat = (cen==False) # (dd['pid']>=1)
			#computes the histogram for each resampling of the file
			for ii, xel in enumerate(X):
				#print ii
				xmin, ymin, zmin, xmax, ymax, zmax = X[ii], Y[ii], Z[ii], X[ii]+Ljk, Y[ii]+Ljk, Z[ii]+Ljk
				sel = (dd['x']>=xmin)&(dd['x']<xmax)&(dd['y']>=ymin)&(dd['y']<ymax)&(dd['z']>=zmin)&(dd['z']<zmax)&(dd[name]>bins[0])&(dd[name]<bins[-1])
				#print len(dd[name][(sel)&(cen)]), len(dd[name][(sel)&(sat)])
				if len(dd[name][(sel)&(cen)])>=1:
					nnC[jj][ii] = n.histogram(dd[name][(sel)&(cen)], bins = bins)[0]
				if  len(dd[name][(sel)&(sat)])>=1:
					nnS[jj][ii] = n.histogram(dd[name][(sel)&(sat)], bins = bins)[0]
			
		f = open(join(output_dir, rootname +"_Central_JKresampling.pkl"),'w')
		cPickle.dump(n.sum(nnC,axis=0),f)
		f.close()
		f = open(join(output_dir,rootname +"_Satellite_JKresampling.pkl"),'w')
		cPickle.dump(n.sum(nnS,axis=0),f)
		f.close()
		n.savetxt(join(output_dir,rootname+"_"+name+"_JKresampling.bins"),n.transpose([bins]))


	def computeSingleDistributionFunctionV2(self, fileList, rootname, name, bins ) :
		"""
		Extracts the distribution of quantity 'name' out of all snapshots of the Multidark simulation.
		Resamples the box in smaller boxes of length Ljk in Mpc/h
		:param ii: index of the snapshot in the list self.snl
		:param name: name of the quantity of interest, mass, velocity.
		:param index: of the quantity of interest in the snapshots.
		:param bins: binning scheme to compute the historgram.
		:param Ljk: length of the resampled box
		:param overlap: allowed overlap between resampled realizations : 1 = no overlap 2 : 50% overlap ... 
		"""		
		output_dir = join(self.wdir,self.boxDir,"properties",name)
		os.system('mkdir '+ output_dir)
		# define boundaries
		nnC = n.zeros((len(fileList),len(bins)-1))
		nnS = n.zeros((len(fileList),len(bins)-1))
		for jj, file in enumerate(fileList):
			#print file
			dd = fits.open(file)[1].data
			cen = (dd['pid']==-1)
			#computes the histogram for each resampling of the file
			sel = (dd['x']>0)&(dd['x']<self.Lbox.value)&(dd['y']>0)&(dd['y']<self.Lbox.value)&(dd['z']>0)&(dd['z']<self.Lbox.value)&(dd[name]>bins[0])&(dd[name]<bins[-1])
			#print len(dd[name][(sel)&(cen)]), len(dd[name][(sel)&(cen==False)])
			if len(dd[name][(sel)&(cen)])>=1:
				nnC[jj] = n.histogram(dd[name][(sel)&(cen)], bins = bins)[0]
			if  len(dd[name][(sel)&(cen==False)])>=1:
				nnS[jj] = n.histogram(dd[name][(sel)&(cen==False)], bins = bins)[0]
			
		f = open(join(output_dir, rootname +"_Central_hist.pkl"),'w')
		cPickle.dump(n.sum(nnC,axis=0),f)
		f.close()
		f = open(join(output_dir,rootname +"_Satellite_hist.pkl"),'w')
		cPickle.dump(n.sum(nnS,axis=0),f)
		f.close()
		n.savetxt(join(output_dir,rootname+"_"+name+".bins"),n.transpose([bins]))

	def computeSingleDistributionFunction(self, ii, name, bins, Mfactor=10. ) :
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
		#print output_dir
		os.system('mkdir '+ output_dir)
		NperBatch = 10000000
		qtyCentral = n.empty(NperBatch)  # 10M array
		qtySat = n.empty(NperBatch)  # 10M array
		#print name, index, output_dir

		fl = fileinput.input(self.snl[ii])
		nameSnapshot = self.snl[ii].split('/')[-1][:-5]
		#print nameSnapshot
		#print name

		countCen,countSat,countFileCen,countFileSat = 0,0,0,0
		
		for line in fl:
			if line[0] == "#" :
				continue

			line = line.split()
			sat_or_cen = float(line[self.columnDict['pid']])
			mv = float(line[self.columnDict['mvir']])
			point = float(line[index])
			if sat_or_cen != -1 and mv > Mfactor * self.Melement and point > 10**bins[0] and point < 10**bins[-1] :
				qtySat[countSat] = point
				countSat+= 1					
				
			if sat_or_cen == -1 and mv > Mfactor * self.Melement and point > 10**bins[0] and point < 10**bins[-1] :
				qtyCentral[countCen] = point
				countCen+= 1					
				
			if countCen == NperBatch :
				#print "cen", qtyCentral
				nnM,bb = n.histogram(n.log10(qtyCentral),bins = bins)
				#print "countCen",countCen
				f = open(join(output_dir, nameSnapshot + "_" + name + "_Central_" + str(countFileCen)+ ".pkl"),'w')
				cPickle.dump(nnM,f)
				f.close()
				countFileCen+= 1
				countCen = 0
				qtyCentral = n.empty(NperBatch)

			if countSat == NperBatch :
				#print "sat", qtySat
				nnM,bb = n.histogram(n.log10(qtySat),bins = bins)
				#print "countSat", countSat
				f = open(join(output_dir, nameSnapshot + "_" + name+ "_Satellite_" + str(countFileSat)+ ".pkl"),'w')
				cPickle.dump(nnM,f)
				f.close()
				countFileSat+= 1
				countSat = 0
				qtySat = n.empty(NperBatch)

		#print "sat2", qtyCentral
		#print "cen2", qtySat
		# and for the last batch :
		nnM,bb = n.histogram(n.log10(qtyCentral[ (qtyCentral>10**bins[0]) & (qtyCentral<10**bins[-1]) ]),bins = bins)
		f = open(join(output_dir, nameSnapshot + "_" + name +"_Central_" + str(countFileCen)+ ".pkl"),'w')
		cPickle.dump(nnM,f)
		f.close()

		nnM,bb = n.histogram(n.log10(qtySat[ (qtySat>10**bins[0]) & (qtySat<10**bins[-1]) ]),bins = bins)
		f = open(join(output_dir, nameSnapshot + "_" + name + "_Satellite_" + str(countFileSat)+ ".pkl"),'w')
		cPickle.dump(nnM,f)
		f.close()
		
		n.savetxt(join(output_dir,name+".bins"),n.transpose([bins]))

	def combinesSingleDistributionFunction(self, ii, name='Vpeak', bins=10**n.arange(0,3.5,0.01), type = "Central" ) :
		"""
		Coombines the outputs of computeSingleDistributionFunction.
		:param ii: index of the snapshot
		:param name: name of the quantity studies
		:param bins: bins the histogram was done with
		:param type: "Central" or "Satellite"
		"""
		output_dir = join(self.wdir,self.boxDir,"properties",name)
		#print output_dir
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
		#print nameA, nameB, indexA, indexB, output_dir

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
				arrA = n.log10(qtyCentral.T[0][(qtyCentral.T[0]>0)])
				arrB = n.log10(qtyCentral.T[1][(qtyCentral.T[1]>0)])
				#print len(arrA), arrA, binsA
				#print len(arrB), arrB, binsB
				nnA,bbA = n.histogram(arrA,bins = binsA)
				nnB,bbB = n.histogram(arrB,bins = binsB)
				dataAB = n.histogram2d(arrA,arrB,bins = [binsA,binsB])
				#print "countCen",countCen
				f = open(join(output_dir, nameSnapshot + "_" + nameA+"-"+nameB + "_Central_" + str(countFileCen)+ ".pkl"),'w')
				cPickle.dump([nnA,nnB,dataAB],f)
				f.close()
				countFileCen+= 1
				countCen = 0
				qtyCentral = n.empty((NperBatch,2))

			if countSat == NperBatch-1 :
				arrAs = n.log10(qtySat.T[0][(qtySat.T[0]>0)])
				arrBs = n.log10(qtySat.T[1][(qtySat.T[1]>0)])
				#print len(arrAs), arrAs, binsA
				#print len(arrBs), arrBs, binsB
				nnA,bbA = n.histogram(arrAs,bins = binsA)
				nnB,bbB = n.histogram(arrBs,bins = binsB)
				dataAB = n.histogram2d(arrAs, arrBs ,bins = [binsA,binsB])
				#print "countSat", countSat
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
		#print "countCen",countCen
		f = open(join(output_dir, nameSnapshot + "_" + nameA+"-"+nameB + "_Central_" + str(countFileCen)+ ".pkl"),'w')
		cPickle.dump([nnA,nnB,dataAB],f)
		f.close()

		nnA,bbA = n.histogram(n.log10(qtySat.T[0]),bins = binsA)
		nnB,bbB = n.histogram(n.log10(qtySat.T[1]),bins = binsB)
		dataAB = n.histogram2d(n.log10(qtySat.T[0]), n.log10(qtySat.T[1]) ,bins = [binsA,binsB])
		#print "countSat", countSat
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




