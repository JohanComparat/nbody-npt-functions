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
# outputs in $BOX_DIR/work_agn or work_cluster
# rewriting agaim the ones with less column numbers

# writes in the catalog dir
#---------------------------
# add stellar masses according to Moster et al. 2013
# to be updated to the Moster et al. 2017 model EMERGE
python $BOX_add_Ms.py
# outputs in $BOX_DIR/work_agn 
# measures the stellar mass function. 
# Is now done in the tabulate duty cycle step
# python measure_SMF.py 

# measure the stellar mass function obtained per snapshot
# and tabulates the duty cyle as a function of stellar mass
# forces the snapshot to reproduce the luminosity function from Bongiorno 2016
python $BOX_tabulate_duty_cycle.py
# outputs in $BOX_DIR/duty_cycle

# add Xray luminosities for AGNs using Bongiorno et al. 2016 and Xray for clusters using Mantz et al. 2016
python $BOX_add_Xray.py
# outputs in $BOX_DIR/work_agn 

#selects active AGNS and write the AGN snapshot in the catalog dir
python $BOX_create_AGN_summary_file.py
# outputs in $BOX_DIR/catalogs/

# add 4MOST targets on top
python MD10_add_4MOST_AGN.py
python MD10_add_4MOST_CLUSTERS_bcg.py
python MD10_add_4MOST_CLUSTERS_members.py
python MD10_add_4MOST_COSMO.py

# 4MOST light cone
# eRosita light cone

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
