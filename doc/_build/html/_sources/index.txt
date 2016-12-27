
Welcome to Skies and Universes python documentation
===================================================

This set of modules helps exploiting the observed and simulated data available through the skies and universes data base http://projects.ift.uam-csic.es/skies-universes/

It currently contains observed galaxy spectra from DEEP2, VVDS and VIPERS. It also hosts parts of the N body simulation: MultiDark and its projections into mock catalogs for emission line galaxies and quasars.
 
Acknowledgment
===========

If you use this code as a ressource to produce a scientific result, please :
 * cite the paper Comparat et al. 2016 https://arxiv.org/abs/1605.02875
 * and acknowledge the skies and universes database as follows : "The skies and universes database is hosted by the UAM/CSIC-IFT cloud and funded by the Spanish MICINNs Consolider-Ingenio 2010 Programme, the MultiDark CSD2009-00064 grant and the MINECO Centro de Excelencia Severo Ochoa Programme under the grants SEV-2012-0249, FPA2012-34694, and the projects AYA2014-60641-C2-1-P and AYA2012-31101."

Install
=====

requirements :
 - python 2.7.8 and its main packages all installable through pip: numpy, scipy, matplotlib, math ...
 - astro dependencies: astropy, pyneb. Both are also installable with pip

git clone https://github.com/JohanComparat/pySU.git

Add all the python folders you can list them like this: ls $PYSU_DIR/*/python) to your pythonpath.

Galaxy surveys
==========

Modules and classes to handle the data
---------------------------------------------

.. toctree::
   :maxdepth: 2 

   GalaxySurveyDEEP2
   GalaxySpectrumDEEP2
   GalaxySurveyVIPERS
   GalaxySpectrumVIPERS
   GalaxySurveyVVDS
   GalaxySpectrumVVDS
   SpectraStacking

Modules and classes to fit models to data
-----------------------------------------------

.. toctree::
   :maxdepth: 2 

   LineFittingLibrary
   LineLuminosityFunction
   ModelSpectraStacks

Support libraries
-----------------------------

.. toctree::
   :maxdepth: 2 

   lib_plot
   filterList
   lineListAir
   MiscellanousFunctionsLibrary

Scripts to construct catalogs
---------------------------------

 * calibrate_DEEP2_spectra : performs the flux calibration of DEEP2 spectra
 * fit_lines_DEEP2_fc_spectra : fits emission lines on the DEEP2 spectra
 * fit_lines_VIPERS_spectra : fits emission lines on the VIPERS spectra
 * fit_lines_VVDSWIDE_spectra : fits emission lines on the VVDS WIDE spectra
 * fit_lines_VVDSDEEP_spectra : fits emission lines on the VVDS DEEP spectra
 * fit_lines_VVDSUDEEP_spectra : fits emission lines on the VVDS UDEEP spectra
 * compute_line_luminosities_DEEP2 : adds line luminosities in a given cosmology
 * compute_line_luminosities_VIPERS : adds line luminosities in a given cosmology + aperture correction
 * compute_line_luminosities_VVDS : idem
 * plotSurveys : produces summary plots related to the surveys

Scripts to measure luminosity functions
----------------------------------------------
In the folder galaxy/bin_LF

 * runLFanalysis : all the steps to reproduce the [OII], Hbeta, [OIII] LF paper, Comparat et al. 2016 https://arxiv.org/abs/1605.02875
 * estimate_line_LF : estimates the emission line luminosity function of the DEEP2, VVDS, and VIPERS surevy

 Scripts to stack spectra
 ---------------------------
 In the folder galaxy/bin_stack

 * stack_spectra_DEEP2 or _VVDSDEEP : stacks the spectra issued from the luminosity function
 
MultiDark N body simulations 
=====================

For exhaustiveinformation on MultiDark, please visit: https://www.cosmosim.org/

Modules and classes 
-------------------------------

.. toctree::
   :maxdepth: 2 

   MultiDark
   LightconeCreation
   HaloSelection
   lib_functions_1pt

Scripts 
-------------
in the folder multidark,
 * bin_DF : scripts related to density field analysis
 * bin_MD : scripts transforming MultiDark database catalogs into light versions to work with
 * bin_onePT : scripts to compute 1-point functions
 * bin_twoPT : scripts to compute 2-point functions
  
  
Stellar population model 
=================
Library in development (caution !), it is the twin of FIREFLY developed in ICG portsmouth and described in Wilkinson et al. 2015 https://arxiv.org/abs/1503.01124 and soon in greater details in Wilkinson et al. in prep.

.. toctree::
   :maxdepth: 2 

   firefly_dust
   firefly_fitter
   firefly_instrument
   firefly_library
   GalaxySpectrumFIREFLY
   StellarPopulationModel

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

