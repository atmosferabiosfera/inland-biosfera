Inland Project's README file.

Package Overview
================

The Integrated Model of Land Surface Processes (Inland) is the land surface package for the Brazilian Earth System Model. The Inland project is coordinated by Marcos Heil Costa (UFV) and Gilvan Sampaio (INPE), with funding from Brazil´s Ministry of Science and Technology programs Rede Clima and INCT.

When completed, Inland will represent a number of land surface processes:
1. Canopy biophysics, including water and energy balance, ice processes, turbulence and aerodynamics
2. Plant physiology processes, including light interception, photosynthesis and respiration, canopy conductance and interactions with soil fertility
3. Soil physics, including energy ans water balance, infiltration, surface runoff and deep drainage 
4. Short-term vegetation dynamics (phenology)
5. Long-term vegetation dynamics, including GPP, NPP, allocation, primary and secondary growth, mortality and competition among plant functional types
6. Distturbances and fires
7. Crop phenology and management, including planting dates, emergence, growth, grain filling, senescence and harvest
8. Soil biogeochemistry, with varying levels of complexity fro carbon, nitrogen and phosphorous
9. Land surface hydrological processes, including rivers, lakes, wetlands and floodplains

The model is forced by meteorological hourly data or regional daily or monthly climate data, or can be coupled to regional or global climate models. 

Inland is based on the IBIS model from SAGE center at the University of Wisconsin-Madison (http://www.sage.wisc.edu/download/IBIS/ibis.html)

For a complete list of people that worked on the INLAND and IBIS projects, see the AUTHORS file. 

Main modifications from IBIS to INLAND
======================================
Released now (2012):
1. All code converted to fortran 90, double precision and dynamic memory allocation
2. Grid and single point versions unified in the same code (different compilations though)
3. Code is multi-thread ready
4. Offline code accepts four frequencies of input data: climate monthly means, monthly time series and daily time series (grid version), hourly time series (single point version)
5. Documentation prepared (the most detailed documentation is available in portuguese only)
6. Uses the GNU Build System

To appear in the next release (expected in 2013):
1. Incorporation of the crop and hydrological codes in the same version (AgroIBIS and THMB)
2. Incorporation of the fire code in the same version
3. Specific pft parameterization for South American biomes
4. Possibility of using Optis, a parameter optimizer for Inland
5. Drivers for coupling to atmospheric models


Main References of the IBIS model
=================================

Foley, J.A. et al., 1996:  An integrated biosphere model of land surface processes, terrestrial carbon balance, and vegetation dynamics.  Global Biogeochemical Cycles, Vol.10, No.4, pages 603-628.

Kucharik, C.J. et al., 2000: Testing the performance of a dynamic global ecosystem model: Water balance, carbon balance and vegetation structure. Global Biogeochemical Cycles 14(3), 795-825.


DIRECTORY TREE OF Inland/IBIS
inland
  |---contrib            Contributions from developers, not required to model execution.
  |---docs               Description of the code.
  |---src                Source code model.
  |---include            Contains parameters for compile-time operations.
  |---data 
     |---offline
       |---grid          Files and parameters to run simulations in GRID mode.
       |  |---input      Boundary conditions and atmospheric forcings.
       |  |---params     Parameter files of canopy, vegetation, soil and fire.
       |  |---conf       Configuration files for the simulation.
       |  |---output     Model output files.
       |
       |---single_point  Files and parameters to run simulations in SINGLE POINT mode.
          |---input      Boundary conditions and atmospheric forcings.
          |---params     Parameter files of canopy, vegetation, soil and fire.
          |---conf       Configuration files for the simulation.
          |---output     Model output files.

More information can be found in the directory 'docs/'

As an effort to keep the model up to date and in the most useful form, we would like to receive your feedback. If you make any substantial changes to the code please let us know. Your change may be incorporated in the next version of the code. The code is under constant development and your suggestion will help us to improve the code. If you experience any problems with the model or have any questions please do not hesitate to contact us at inland@inpe.br

August 27th, 2012.

Regards,
Dr. Gilvan Sampaio de Oliveira - gilvan.sampaio@inpe.br (http://ccst.inpe.br)
Dr. Marcos Heil Costa - mhcosta@ufv.br (http://madeira.dea.ufv.br/)
and the Inland R&D team - inland@inpe.br
