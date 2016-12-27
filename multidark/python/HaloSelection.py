
"""
Library of function to create halo catalogs matched to a density.

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

"""
import random
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as p
import astropy.io.fits as fits
import numpy as n
from scipy.interpolate import interp1d
import scipy.stats as st
import os
from os.path import join

class MultiDarkMock:
    """
    :param hdu: hdu of the lightcone
    :param area: area in deg2
    :param mockOutput_dir: directory where to output mocks
    :param mockName: name of the file where to save the mock
    :param zmin: minimum redshift array defining the bins
    :param zmax: maximum redshift array defining the bins
    :param nGal_Deg2: number of galaxies per square degrees 
    """
    def __init__(self,hdu, area, mockOutput_dir, mockName, zmin, zmax, nGal_Deg2 ):
        self.hdu = hdu
        self.area = area
        self.mockOutput_dir = mockOutput_dir
        self.mockName = mockName
        self.zmin = zmin
        self.zmax = zmax
        self.nGal_Deg2 = nGal_Deg2
        
    def initialize(self):
        """
        Initializes the procedure by putting into memroy arrays of central and satellites.
        """
        # derived numbers common to all SHAMs.
        self.cen = (self.hdu[1].data['pid'] ==  -1)
        self.sat = (self.cen ==  False)
        self.Nhalos = self.hdu[1].header['NAXIS2']
        self.IDh = n.arange(self.Nhalos)
        self.nGal = n.array([ int(el * self.area)  for el in self.nGal_Deg2 ])
        #self.nGal_to_z = interp1d(nGal_Deg2,(self.zmax+self.zmin)/2.)
        #function to slice by redshift 
        self.slice_Z = lambda z1, z2 : (self.hdu[1].data['z_redshift_space'] >= z1) & ( self.hdu[1].data['z_redshift_space'] < z2)

    def write_catalog_ascii(self):
        """Writes the obtained mock catalog for clustering estimation: just ra, dec and redshift. """
        print "writes ascii catalog :", self.mockName
        outPutFileName = join( self.mockOutput_dir, self.mockName + "_radecz.cat" )
        self.raMock = self.hdu[1].data['ra'][self.idSel]
        self.decMock = self.hdu[1].data['dec'][self.idSel]
        self.zMock = self.hdu[1].data['z_redshift_space'][self.idSel]
        n.savetxt(outPutFileName, n.transpose([ self.raMock, self.decMock, self.zMock]),fmt = '%.8f %.8f %.5f')

    def write_full_catalog_fits(self):
        """Writes the obtained with all the columns from the parent lightcone catalog."""
        print "writes fits catalog :", self.mockName
        tbhdu = fits.BinTableHDU.from_columns( self.hdu[1].columns )
        tbhdu.data = tbhdu.data[self.idSel]
        prihdu = fits.PrimaryHDU(header = self.hdu[0].header)
        thdulist = fits.HDUList([prihdu, tbhdu])
        outPutFileName = join(self.mockOutput_dir,self.mockName+"_allCols.fits")
        os.system('rm -rf '+ outPutFileName)
        thdulist.writeto(outPutFileName)

    def get_distrib_QTY(self, colN, z1, z2):
        """Computes the cumulative histogram of a column for halos in the range z1, z2.
        :param colN: name of the column you want to take the histogram.
        :param z1: minimum redshift
        :param z2: maximum redshift
        """
        zsel = self.slice_Z( z1, z2)
        IDhz = self.IDh[zsel] # all ids in this redshift bin
        QTY = self.hdu[1].data[colN][zsel] # all QTY in this redshift bin
        nn,bb,pp = p.hist(QTY,cumulative = True,bins = len(QTY)/100)
        p.clf()
        print len(IDhz), "halos with ",z1, "<z<", z2
        return IDhz,QTY,nn,bb

    def select_sham(self, nGal_perbin, IDhz, QTY, nn, bb):
        """
        Returns the ids corresponding to a given density.
        :param nGal_perbin: number of galaxies to be selected
        :param IDhz: parent ID distribution
        :param QTY: quantity to select halos on
        :param nn: cumulative distribution o QTY
        :param bb: bins of the cumulative distribution
        """
        mfc = interp1d(nn, (bb[:-1]+bb[1:])/2.)
        QTYmax = mfc(len(QTY))
        QTYmin = mfc(len(QTY)-nGal_perbin)
        qsel = (QTY>QTYmin)&(QTY<= QTYmax)
        IDhzq = IDhz[qsel]
        print "N to be selected:",nGal_perbin,", Nselected:",len(IDhzq)
        return IDhzq

    def make_sham_catalog(self, colN='mvir'):
        """
        Creates lists of ids of halos corresponding to the density given in the n(z). 
        For every bin of redshift, it gets the distribution o fthe column of interest and matches to the density of galaxies in the NZ given.
        Then provides a column of ids extracted from teh lightcone.
        :param colN: name of the column you wish to work on for the sham.
        """
        ids = []
        for ii in range(len(self.zmin)): 
            print "gets all halos for ", self.zmin[ii], "<z<", self.zmax[ii], "with col5 to mock ", self.nGal[ii], " galaxies." 
            IDhz, QTY, nn, bb = self.get_distrib_QTY( colN, self.zmin[ii], self.zmax[ii])
            ids.append( self.select_sham(self.nGal[ii],IDhz, QTY, nn,bb)) 

        self.idSel = n.hstack(( ids ))
        self.NhaloMock = len((self.idSel).nonzero()[0])

    def select_shamIncomplete(self, incompFactor, nGal_perbin, IDhz, QTY, nn, bb):
        """
        Returns the ids corresponding to a given density and an incompleteness factor.
        :param nGal_perbin: number of galaxies to be selected
        :param IDhz: parent ID distribution
        :param QTY: quantity to select halos on
        :param nn: cumulative distribution o QTY
        :param bb: bins of the cumulative distribution
        :param incompFactor: incompleteness factor compared to the max of QTY : max(QTY)/incompFactor will be set as the max of the extracted distribution.
        """
        mfc = interp1d(nn,(bb[1:]+bb[:-1])/2.)
        mfcInv = interp1d((bb[1:]+bb[:-1])/2.,nn)
        QTYmaxAll = mfc(len(QTY))/incompFactor
        Nmax = mfcInv(QTYmaxAll)
        QTYmax = mfc(Nmax)
        QTYmin = mfc(Nmax-nGal_perbin)
        qsel = (QTY>QTYmin)&(QTY<= QTYmax)
        IDhzq = IDhz[qsel]
        return IDhzq
        
    def make_shamIncomplete_catalog(self, colN, incompletenessFactor ):
        """
        Creates lists of ids of halos corresponding to the density given in the n(z). 
        For every bin of redshift, it gets the distribution o fthe column of interest and matches to the density of galaxies in the NZ given.
        Then provides a column of ids extracted from teh lightcone.
        :param colN: name of the column you wish to work on for the sham.
        """
        ids = []
        for ii in range(len(self.zmin)):
            print "gets all halos for ", self.zmin[ii], "<z<", self.zmax[ii], "with col5 to mock ", self.nGal[ii], " galaxies." 
            IDhz, QTY, nn, bb = get_distrib_QTY( hdu, colN, self.zmin[ii], self.zmax[ii] )
            ids.append( self.select_shamIncomplete( incompletenessFactor[ii], self.nGal[ii], IDhz, QTY, nn, bb ) ) 

        self.idSel = n.hstack(( ids ))
        self.NhaloMock = len((self.idSel).nonzero()[0])

    def select_shamMAX(self,QTY_max, nGal_perbin,IDhz, QTY, nn,bb):
        """
        Returns the ids corresponding to a given density with a maximum in the QTY.
        :param nGal_perbin: number of galaxies to be selected
        :param IDhz: parent ID distribution
        :param QTY: quantity to select halos on
        :param nn: cumulative distribution o QTY
        :param bb: bins of the cumulative distribution
        :param QTY_max: the max of QTY is set to QTY_max.
        """
        mfc = interp1d(nn,(bb[1:]+bb[:-1])/2.)
        mfcInv = interp1d((bb[1:]+bb[:-1])/2.,nn)
        Nmax = mfcInv(QTY_max)
        QTYmax = mfc(Nmax)
        QTYmin = mfc(Nmax-nGal_perbin)
        qsel = (QTY>QTYmin)&(QTY<= QTYmax)
        IDhzq = IDhz[qsel]
        return IDhzq

    def make_shamMAX_catalog(self, colN, maxQTY ):
        """
        Creates lists of ids of centrals and satellite galaxies corresponding to the density given in the n(z). 
        For every bin of redshift, it gets the distribution o fthe column of interest and matches to the density of galaxies in the NZ given.
        Then provides a column of ids extracted from teh lightcone.
        :param colN: name of the column you wish to work on for the sham.
        """
        ids = []
        for ii in range(len(self.zmin)):
            print "gets all halos for ", self.zmin[ii], "<z<", self.zmax[ii], "with col5 to mock ", self.nGal[ii], " galaxies." 
            IDhz, QTY, nn, bb = get_distrib_QTY( hdu, colN, self.zmin[ii], self.zmax[ii] )
            ids.append( self.select_shamMAX( maxQTY[ii], self.nGal[ii], IDhz, QTY, nn, bb ) ) 

        self.idSel = n.hstack(( ids ))
        self.NhaloMock = len((self.idSel).nonzero()[0])

    def select_Gaussian(self, meanQTY, scatterQTY, nGal_perbin, IDhz, QTY):
        """
        Creates lists of ids of centrals and satellite galaxies corresponding to the density given in the n(z) and to a gaussian distribution . 
        For every bin of redshift, it gets the distribution o fthe column of interest and matches to the density of galaxies in the NZ given.
        Then provides a column of ids extracted from teh lightcone.
        :param colN: name of the column you wish to work on for the sham.
        :param meanQTY: mean of the distribution
        :param scatterQTY: scatter of the distribution
        :param nGal_perbin: total number of galaxies in this bins to mock
        :param IDhz: IDs of the halos in this bin
        :param QTY: array of the column to do the match on, mass, velocity, ...
        """
        # constructs the QTY intervals around the distribution
        expected_cdf = lambda x : st.norm.cdf(x, loc = meanQTY, scale = scatterQTY)
        interval  =  [ meanQTY - 9 * scatterQTY , meanQTY + 9 * scatterQTY]
        xs = n.arange(interval[0],interval[1],(interval[1]-interval[0])/1000.)
        out = expected_cdf(xs)
        expected_cdf_inv = interp1d(out,xs)
        boundaries = n.hstack((expected_cdf_inv(0.01),expected_cdf_inv(n.arange(0.1,0.91,0.1)), interval[1]))
        # gets the number of halos to select
        expected_cdf_tot = lambda x : nGal_perbin * st.norm.cdf(x, loc = meanQTY, scale = scatterQTY)
        Up = expected_cdf_tot(boundaries[1:])
        Low = n.hstack(( 0., expected_cdf_tot(boundaries[1:])[:-1] ))
        N2select = Up-Low
        print N2select,Up,Low
        # select in mass in the box
        qsels = n.array([ (QTY>boundaries[ii])&(QTY<= boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
        IDhzqAll = n.array([ IDhz[qs] for qs in qsels ])
        # random downsample to the N2select in each bin
        i = 0
        ids_selected = []
        for arr in IDhzqAll:
            random.shuffle(arr)
            ids_selected.append(arr[:N2select[i]])
            i+= 1

        ids_selected = n.hstack(( n.array(ids_selected) ))
        return ids_selected

    def make_Gaussian_catalog(self, colN, means, scatters):
        """
        Creates lists of ids of centrals and satellite galaxies corresponding to the density given in the n(z).
        :param colN: name of the column you construct the catalog with
        :param means: means of the Gaussians, array the same length of the redshift bin
        :param scatters: scatters of the Gaussians, array the same length of the redshift bin
        """
        ids = []
        for ii in range(len(self.zmin)):
            print "gets all halos for ", self.zmin[ii], "<z<", self.zmax[ii], "with col5 to mock ", self.nGal[ii], " galaxies." 
            IDhz, QTY, nn, bb = get_distrib_QTY( hdu, colN, self.zmin[ii], self.zmax[ii] )
            ids.append( self.select_Gaussian( means[ii], scatters[ii], self.nGal[ii], IDhz, QTY ) ) 

        self.idSel = n.hstack(( ids ))
        self.NhaloMock = len((self.idSel).nonzero()[0])
        
    def get_distrib_QTY_cen(self, colN, z1, z2):
        """Computes the cumulative histogram of a column for central halos in the range z1, z2.
        :param colN: name of the column you want to take the histogram.
        :param z1: minimum redshift
        :param z2: maximum redshift
        """
        zsel = self.slice_Z(z1, z2) & (self.cen)
        IDhz = self.IDh[zsel] # all ids in this redshift bin
        QTY = self.hdu[1].data[colN][zsel] # all QTY in this redshift bin
        nn,bb,pp = p.hist(QTY,cumulative = True,bins = len(QTY)/100)
        p.clf()
        return IDhz,QTY,nn,bb

    def get_distrib_QTY_sat(self, colN, z1, z2):
        """Computes the cumulative histogram of a column for satellite halos in the range z1, z2.
        :param colN: name of the column you want to take the histogram.
        :param z1: minimum redshift
        :param z2: maximum redshift
        """
        zsel = self.slice_Z(z1, z2) & (self.sat)
        IDhz = self.IDh[zsel] # all ids in this redshift bin
        QTY = self.hdu[1].data[colN][zsel] # all QTY in this redshift bin
        nn,bb,pp = p.hist(QTY,cumulative = True,bins = len(QTY)/100)
        p.clf()
        return IDhz,QTY,nn,bb

    def select_GaussianFsat(self,meanQTY,scatterQTY,fsat, nGal_perbin, IDhz_c, QTY_c, IDhz_s, QTY_s ):
        """
        Extracts the ids of halos to create a mock with a Gaussian distribution.
        :param colN: name of the column you wish to work on for the sham.
        :param meanQTY: mean of the distribution
        :param scatterQTY: scatter of the distribution
        :param fsat: fraction of satellite in this bin
        :param nGal_perbin: total number of galaxies in this bins to mock
        :param IDhz_c: IDs of the central halos in this bin
        :param QTY_c: column to do the match on, mass, velocity, ... for the central halos
        :param IDhz_s: IDs of the satellite halos in this bin
        :param QTY_s: column to do the match on, mass, velocity, ... for the satellite halos
        """
        nSat = int(nGal_perbin*fsat)
        print "satellites",nGal_perbin,nSat,fsat,meanQTY,scatterQTY
        # constructs the QTY intervals around the distribution 
        expected_cdf = lambda x : st.norm.cdf(x, loc = meanQTY, scale = scatterQTY)
        interval  =  [ meanQTY - 9 * scatterQTY , meanQTY + 9 * scatterQTY]
        xs = n.arange(interval[0],interval[1],(interval[1]-interval[0])/1000.)
        out = expected_cdf(xs)
        expected_cdf_inv = interp1d(out,xs)
        boundaries = n.hstack((expected_cdf_inv(0.01),expected_cdf_inv(n.arange(0.1,0.91,0.1)), interval[1]))
        # gets the number of halos to select the SAT
        expected_cdf_s  =  lambda x : nSat * st.norm.cdf(x, loc = meanQTY, scale = scatterQTY)
        Up_s  =  expected_cdf_s(boundaries[1:])
        Low_s  =  n.hstack(( 0., expected_cdf_s(boundaries[1:])[:-1] ))
        N2select_s  =  Up_s-Low_s
        # select in mass in the box
        qsels_s = n.array([ (QTY_s>boundaries[ii])&(QTY_s<= boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
        IDhzqAll_s = n.array([ IDhz_s[qs] for qs in qsels_s ])

        # random downsample to the N2select in each bin
        i = 0
        ids_selected_s = []
        for arr2 in IDhzqAll_s:
            random.shuffle(arr2)
            #print len(arr2),int(N2select_s[i])
            ids_selected_s.append(arr2[:int(N2select_s[i])])
            i+= 1

        id_s = n.hstack((n.array(ids_selected_s)))

        nSatReal = len(id_s)
        nCen = nGal_perbin-nSatReal
        print "centrals", nGal_perbin,nSat,nCen,fsat,meanQTY,scatterQTY

        # gets the number of halos to select the CEN, compatible with the sat fraction to get the right density.
        print "centrals"
        expected_cdf_c  =  lambda x : nCen * st.norm.cdf(x, loc = meanQTY, scale = scatterQTY)
        Up_c  =  expected_cdf_c(boundaries[1:])
        Low_c  =  n.hstack(( 0., expected_cdf_c(boundaries[1:])[:-1] ))
        N2select_c  =  Up_c-Low_c
        # select in mass in the box
        qsels_c = n.array([ (QTY_c>boundaries[ii])&(QTY_c<= boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
        IDhzqAll_c = n.array([ IDhz_c[qs] for qs in qsels_c ])

        # random downsample to the N2select in each bin
        i = 0
        ids_selected_c = []
        for arr in IDhzqAll_c:
            random.shuffle(arr)
            #print len(arr),int(N2select_c[i])
            ids_selected_c.append(arr[:int(N2select_c[i])])
            i+= 1

        id_c = n.hstack((n.array(ids_selected_c)))

        ids_selected = n.hstack((id_c,id_s ))
        print len(id_c),len(id_s),len(ids_selected)
        return ids_selected


    def make_GaussianFsat_catalog(self, colN, means, scatters, fsats):
        """
        Creates lists of ids of centrals and satellite galaxies corresponding to the density given in the n(z).
        :param colN: name of the column you construct the catalog with
        :param means: means of the Gaussians, array the same length of the redshift bin
        :param scatters: scatters of the Gaussians, array the same length of the redshift bin
        :param fsats: fractions of satellite, array the same length of the redshift bin
        """
        ids = []
        for ii in range(len(self.zmin)):
            print "gets all halos for ",self.zmin[ii],"<z<",self.zmax[ii], "with col5 to mock ", self.nGal[ii], " galaxies." 
            IDhz_c,QTY_c,nn_c,bb_c = self.get_distrib_QTY_cen( colN, z1=self.zmin[ii], z2=self.zmax[ii])
            IDhz_s,QTY_s,nn_s,bb_s = self.get_distrib_QTY_sat( colN, z1=self.zmin[ii], z2=self.zmax[ii])
            ids.append( self.select_GaussianFsat( means[ii], scatters[ii], fsats[ii], self.nGal[ii], IDhz_c, QTY_c, IDhz_s, QTY_s  ) ) 

        self.idSel = n.hstack(( ids ))
        self.NhaloMock = len((self.idSel).nonzero()[0])

    def select_LogNorm(self, meanQTY, scatterQTY, nGal_perbin,IDhz, QTY, nn,bb):
        """
        Creates lists of ids of centrals and satellite galaxies corresponding to the density given in the n(z) and to a gaussian distribution . 
        For every bin of redshift, it gets the distribution o fthe column of interest and matches to the density of galaxies in the NZ given.
        Then provides a column of ids extracted from teh lightcone.
        :param colN: name of the column you wish to work on for the sham.
        :param meanQTY: mean of the distribution
        :param scatterQTY: scatter of the distribution
        :param nGal_perbin: total number of galaxies in this bins to mock
        :param IDhz: IDs of the halos in this bin
        :param QTY: array of the column to do the match on, mass, velocity, ...
        """
        # constructs the QTY intervals around the distribution
        expected_cdf = lambda x : st.lognorm.cdf(x, meanQTY, scatterQTY)
        interval  =  [ meanQTY - 9 * scatterQTY , meanQTY + 9 * scatterQTY]
        xs = n.arange(interval[0],interval[1],(interval[1]-interval[0])/1000.)
        out = expected_cdf(xs)
        expected_cdf_inv = interp1d(out,xs)
        boundaries = n.hstack((expected_cdf_inv(0.01),expected_cdf_inv(n.arange(0.1,0.91,0.1)), interval[1]))
        # gets the number of halos to select
        expected_cdf_tot = lambda x : nGal_perbin * st.lognorm.cdf(x, meanQTY, scatterQTY)
        Up = expected_cdf_tot(boundaries[1:])
        Low = n.hstack(( 0., expected_cdf_tot(boundaries[1:])[:-1] ))
        N2select = Up-Low
        #print N2select,Up,Low
        # select in mass in the box
        qsels = n.array([ (QTY>boundaries[ii])&(QTY<= boundaries[ii+1]) for ii in range(len(boundaries)-1) ])
        IDhzqAll = n.array([ IDhz[qs] for qs in qsels ])
        # random downsample to the N2select in each bin
        i = 0
        ids_selected = []
        for arr in IDhzqAll:
            random.shuffle(arr)
            ids_selected.append(arr[:N2select[i]])
            i+= 1

        ids_selected = n.hstack(( n.array(ids_selected) ))
        return ids_selected

    def make_LogNorm_catalog(self, colN, means, scatters):
        """
        Creates lists of ids of centrals and satellite galaxies corresponding to the density given in the n(z).
        :param colN: name of the column you construct the catalog with
        :param means: means of the Gaussians, array the same length of the redshift bin
        :param scatters: scatters of the Gaussians, array the same length of the redshift bin
        """
        ids = []
        for ii in range(len(self.zmin)):
            print "gets all halos for ", self.zmin[ii], "<z<", self.zmax[ii], "with col5 to mock ", self.nGal[ii], " galaxies." 
            IDhz, QTY, nn, bb = get_distrib_QTY( hdu, colN, self.zmin[ii], self.zmax[ii] )
            ids.append( self.select_LogNorm( means[ii], scatters[ii], self.nGal[ii], IDhz, QTY ) ) 

        self.idSel = n.hstack(( ids ))
        self.NhaloMock = len((self.idSel).nonzero()[0])

    def create_random_catalog(self, factor = 5., dz=0.025 ):
        """Writes a random catalog"""
        self.nRandom = int(self.NhaloMock * factor )
        raR = n.random.uniform(n.min(self.raMock), n.max(self.raMock), self.nRandom )
        decR = n.random.uniform(n.min(self.decMock), n.max(self.decMock), self.nRandom )
        z1=n.arange(n.min(self.zMock)-0.1, n.max(self.zMock)+0.1, dz)
        nn,bb,pp=p.hist(self.zMock, bins=z1)
        nz=interp1d((z1[1:]+z1[:-1])/2.,factor*nn)
        zs=n.arange(n.min(self.zMock), n.max(self.zMock), dz)
        rdsz=[]
        for i in range(len(zs)-1):
            inter=n.random.uniform(low=zs[i], high=zs[i+1], size=int(2* nz( zs[i]+dz/2. )))
            rdsz.append(inter)

        rds=n.hstack((rdsz))
        n.random.shuffle(rds)
        selRDS=(n.random.rand(len(raR))<float(self.nRandom)/len(raR))
        RR=rds[:len(raR[selRDS])]
        print "N final",len(raR[selRDS])
        outPutFileName = join( self.mockOutput_dir, self.mockName + "_random.cat" )
        n.savetxt(outPutFileName,n.transpose([raR[selRDS],decR[selRDS],RR]),fmt='%.8f %.8f %.5f')
        raR,decR,RR=0,0,0
        
    def writeClusteringParamFile(self,type,decade=""):
        """ Writes the clustering commands that command the CUTE code, see Alonso et al. 2012 https://arxiv.org/abs/1210.1833
        :param type: monopole or angular or ...
        :param decade: string suffix that is appended if you study different scales (decades) _d1, _d2, _d3 are used for the angular clustering."""
        f=open(join( self.mockOutput_dir, self.mockName +".param2PCF_"+type+decade),'a')
        f.write("data_filename= "+join( self.mockOutput_dir, self.mockName + "_radecz.cat" )+" \n")
        f.write("random_filename= "+join( self.mockOutput_dir, self.mockName + "_random.cat" )+" \n")
        f.write("input_format= 2 \n")
        f.write("mask_filename= 'none' \n")
        f.write("z_dist_filename= 'none' \n")
        f.write("output_filename= "+join( self.mockOutput_dir, self.mockName )+"_2PCF_"+type+decade+".dat \n")
        f.write("num_lines= all \n")
        f.write("corr_type= "+type+" \n")
        f.write("corr_estimator= LS \n")
        f.write("np_rand_fact= 5 \n")
        f.write("omega_M= 0.307115 \n")
        f.write("omega_L= 0.692885 \n")
        f.write("w= -1 \n")
        f.write("radial_aperture= 1 \n")
        f.write("use_pm= 0 \n")
        f.write("n_pix_sph= 2048 \n")
        f.close()

    def compute_clustering(self):
        """ Runs the CUTE code to estimate clustering using the LS estimator. """
        os.system("/home2/jcomparat/code/CUTE-1.1A1/CUTE/CUTE "+join( self.mockOutput_dir, self.mockName +".param2PCF_angular_d1"))
        os.system("/home2/jcomparat/code/CUTE-1.1A2/CUTE/CUTE "+join( self.mockOutput_dir, self.mockName +".param2PCF_angular_d2"))
        os.system("/home2/jcomparat/code/CUTE-1.1A3/CUTE/CUTE "+join( self.mockOutput_dir, self.mockName +".param2PCF_angular_d3"))
        os.system("/home2/jcomparat/code/CUTE-1.1M/CUTE/CUTE "+join( self.mockOutput_dir, self.mockName +".param2PCF_monopole"))

    def compare_clustering_data_mock(self, w_data, xi_data, theta_min_chi2 = -2.3,  theta_max_chi2= -1.5, w_bins=15, s_min_chi2=0.5, s_max_chi2=1.2, s_bins=10):
        """Compares the clustering of the mock catalog and the clustering of the data.
        :param w_data: angular clustering from the data [x, y, yErr].
        :param xi_data: monopole clustering from the data [x, y, yErr]. """
        ths = n.logspace(theta_min_chi2,theta_max_chi2,w_bins)
        ss = n.logspace(s_min_chi2,s_max_chi2,s_bins)
        
        # loads the angular clustering from the mock
        xx0, yy0, y2E = n.loadtxt( join( self.mockOutput_dir, self.mockName )+"_2PCF_"+"angular"+"_d3"+".dat",unpack=True,usecols = (0,1,2))
        xx1, yy1, y1E= n.loadtxt( join( self.mockOutput_dir, self.mockName )+"_2PCF_"+"angular"+"_d2"+".dat",unpack=True,usecols = (0,1,2))
        w_M = interp1d( n.hstack((xx0,xx1[1:])), n.hstack((yy0,yy1[1:])) ) #,yy2[1:]))

        # loads the monopole from the mock
        s_M_a, xi_M_a = n.loadtxt( join( self.mockOutput_dir, self.mockName )+"_2PCF_"+"monopole"+".dat",unpack=True,usecols = (0,1))
        xi_M = interp1d( s_M_a, xi_M_a )
        # loads the monopole from the mock

        #s_selection_data=( xi_data[0] > s_min_chi2 ) & ( xi_data[0] < s_max_chi2 ) & (xi_data[1] > 2 * xi_data[2])
        #theta_selection_data = ( w_data[0] > theta_min_chi2 ) & ( w_data[0] < theta_max_chi2 ) & (w_data[1] > 2 * w_data[2])

        xi_D = interp1d( xi_data[0], xi_data[1]) #[s_selection_data], xi_data[1][s_selection_data] )
        xi_D_err = interp1d( xi_data[0], xi_data[2]) #[s_selection_data], xi_data[2][s_selection_data] )

        w_D = interp1d( w_data[0], w_data[1]) #[theta_selection_data], w_data[1][theta_selection_data] )
        w_D_err = interp1d( w_data[0], w_data[2]) #[theta_selection_data], w_data[2][theta_selection_data] )

        chi2Wr = n.sum((w_D(ths) - w_M(ths))**2. / w_D_err(ths)**2) /len(ths)
        chi2Xr = n.sum((xi_D(ss) - xi_M(ss))**2. / xi_D_err(ss)**2) /len(ss)
        
        return chi2Wr, chi2Xr
        