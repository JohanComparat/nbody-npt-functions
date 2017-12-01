# verification of the implementation
------------------------------------

# for AGNs
# reproduce the plots from Bongiorno et al. 2010
# writes plots in data/eRoMok/
python reproduce_bongiorno_2010.py


BOX = MD10 or NUGC

#writes a summary file containing all information for each snapshot
#---------------------------
python $BOX_box-get-header.py

# rewrite rockstar ascii catalogs in smaller fits files with 20e6 lines each + halo mass cut.
python $BOX-write-clusterFiles.py # for the cluster calculations
python $BOX-write-smallFile.py # for the AGN calculations
# outputs in /work_agn or work_cluster
python2.7 MD10-check-small-file-1pt-fun.py
python2.7 MD10-check-small-file-1pt-fun-plots.py
# outputs in wwwDir/eRoMok/
# ALL look OK

# writes in the catalog dir
#---------------------------
# add stellar masses according to Moster et al. 2013
# to be updated to the Moster et al. 2017 model EMERGE
python2.7 MD10_add_Ms.py # (and all other scripts Ms_?) 
python2.7 MD10-check-MS-file-1pt-fun.py 
python2.7 MD10-check-MS-file-1pt-fun-plots.py   
python2.7 plot_SMHMR.py

# measures the stellar mass function. 
# Is now done in the tabulate duty cycle step 
python2.7 measure_SMF.py 
# outputs in $BOX_DIR/duty_cycle
python2.7 plot_SMF.py
# outputs in os.path.join(os.environ['MD10'],"results","stellar_mass_function", "images")

#########################################33
# tabulates the duty cyle as a function of stellar mass
# forces the snapshot to reproduce the luminosity function from Bongiorno 2016
python2.7 MD10_tabulate_duty_cycle.py
# outputs in $BOX_DIR/duty_cycle

# add a flag for the activity to reproduce the host galaxy stellar mass function
# output in $MD10/work_agn/*_DC.fits
python2.7 MD10_add_AGN_activity.py

# add LSAR for AGNs using Bongiorno et al. 2016 
# output in $MD10/work_agn/*_LSAR.fits
python2.7 MD10_add_LSAR.py
python2.7 MD10_add_LSAR_1.py
python2.7 MD10_add_LSAR_2.py
python2.7 MD10_add_LSAR_3.py
python2.7 MD10_add_LSAR_4.py

# create a single AGN file per snapshot
# output in $MD10/catalogs/*_LSAR.fits
python2.7 MD10_create_AGN_summary_file.py

# add obscuration following Buchner et al. 2016
# output in $MD10/work_agn/*_NH.fits
python2.7 MD10_add_AGN_obscuration.py

# add Xray luminosities for AGNs using Bongiorno et al. 2016 and Xray for clusters using Mantz et al. 2016
python MD10_add_Xray.py
python MD10_add_Xray_1.py
python MD10_add_Xray_2.py
python MD10_add_Xray_3.py
python MD10_add_Xray_4.py
python MD10_add_Xray_5.py

# outputs in $BOX_DIR/work_agn 

#selects active AGNS and write the AGN snapshot in the catalog dir
python MD10_create_AGN_summary_file.py
# outputs in $BOX_DIR/catalogs/

# add 4MOST targets on top
python MD10_add_4MOST_AGN.py
python MD10_add_4MOST_CLUSTERS_bcg.py
python MD10_add_4MOST_CLUSTERS_members.py
python MD10_add_4MOST_COSMO.py

# 4MOST light cone
# eRosita light cone
python measure_SMF.py
python plot_SMF.py

python measure_HMF_tracer.py
python plot_HMF.py

python MD10-pie-plot.py

# AFTER LC is created
python MD10_select_AGN_lightcone_eRositaFlux_per_file.py
# read and writes here $MD10/light-cone/

# TB UPDATE FROM HERE ON
# plots and results 
python plot_AGN_HGMF_duty_cycle.py
python plot_slice_simulation.py

#-------------------------------------------------------------------------
python plot_cluster_scaling_relations.py
python plot_LFX.py
python plot_Lxcut_Mhalo.py
python plot-Ms-xray.py

python test.py
