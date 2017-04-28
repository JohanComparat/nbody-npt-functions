# verification of the implementation
------------------------------------

# for AGNs
# reproduce the plots from Bongiorno et al. 2010
# writes plots in data/eRoMok/
python reproduce_bongiorno_2010.py


#writes in the snapshot dir
---------------------------

# rewrite catalogs in fits files
python MD04-write-smallFile.py
python MD10-write-smallFile.py
python MD25-write-smallFile.py


# writes in the catalog dir
---------------------------
# this code it to be double checked and merged into a class ...

# add stellar masses according to Moster et al. 2013
python add_Ms.py

# measure the stellar mass function obtained
python measure_SMF.py

# tabulate the duty cyle as a function of stellar mass
python tabulate_duty_cycle.py
python plot_AGN_HGMF_duty_cycle.py

# add Xray luminosities for AGNs
python add_Xray.py


#-------------------------------------------------------------------------
python plot_cluster_scaling_relations.py
python plot_LFX.py
python plot_Lxcut_Mhalo.py
python plot-Ms-xray.py

python test.py
