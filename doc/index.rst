Welcome to Nbody n-point function documentation
===============================================

This set of modules helps exploiting the simulated data available through the skies and universes data base http://projects.ift.uam-csic.es/skies-universes/

It currently hosts parts of the N body simulation: MultiDark and its projections into mock catalogs for emission line galaxies and quasars.
 
Acknowledgment
==============

If you use this code as a ressource to produce a scientific result, please :
 * cite the paper Comparat et al. 2017 https://arxiv.org/abs/1702.01628
 * and acknowledge the skies and universes database as follows : "The skies and universes database is hosted by the UAM/CSIC-IFT cloud and funded by the Spanish MICINNs Consolider-Ingenio 2010 Programme, the MultiDark CSD2009-00064 grant and the MINECO Centro de Excelencia Severo Ochoa Programme under the grants SEV-2012-0249, FPA2012-34694, and the projects AYA2014-60641-C2-1-P and AYA2012-31101."

Install
=======

requirements :
 - python 2.7.8 and its main packages all installable through pip: numpy, scipy, matplotlib, math ...
 - astro dependencies: astropy, pyneb, icrar hmf package. All are also installable with pip

git clone https://github.com/JohanComparat/nbody-npt-functions.git

Add all the python folders you can list them like this: ls $NBODY_NPT_DIR/*/python) to your pythonpath.

MultiDark N body simulations 
============================

For exhaustive information on MultiDark, please visit: https://www.cosmosim.org/

Modules and classes 
-------------------

.. toctree::
   :maxdepth: 2 

   MultiDark
   LightconeCreation
   HaloSelection
   lib_functions_1pt
   XrayLuminosity
   ClusterScalingRelations
   StellarMass


Scripts 
-------

in the folder multidark,
 * bin_DF : scripts related to density field analysis
 * bin_MD : scripts transforming MultiDark database catalogs into lighter versions to work with
 * bin_onePT : scripts to compute 1-point functions
 * bin_twoPT : scripts to compute 2-point functions
  
  
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

