# FLUXOS-OVERLAND
## Table of Contents
* [Introduction](#introduction)
* [Branches](#branches)
* [Compiling](#compiling)
* [Execution](#execution)
* [Visualization](#visualization)
* [Supporting Scripts](#supporting-scripts)
* [Working Example](#working-example)

## Introduction
* Soure code for the FLUXOS-OVERLAND model. The original code (named as FLUXOS) was written in Fortran and consisted of the coupling of 2dmb, +QeS2, MODFLOW and MT3DMS.

* Modifications (from FLUXOS to FLUXOS-OVERLAND):
	* Converted to C++
	* Uses: Armadillo template-based C++ library for linear algebra 
	* Removed: MODFLOW and MT3DMS (currently there is no baseflow)
	* WINTRA algorithm was integrated for calculation of runoff-soil interactions and nutrient release: https://onlinelibrary.wiley.com/doi/full/10.1002/hyp.11346

* Reading material:
	* Theoretical background (original FLUXOS):
		* EMS paper: https://www.sciencedirect.com/science/article/pii/S1364815216306193?via%3Dihub
		* PhD thesis: https://scholarbank.nus.edu.sg/handle/10635/124183
	* Applications (original FLUXOS):
		* STC paper: https://www.sciencedirect.com/science/article/pii/S0169772216300948?via%3Dihub
		* JCH paper: https://www.sciencedirect.com/science/article/pii/S0169772216300948?via%3Dihub
		* JAWRA: https://onlinelibrary.wiley.com/doi/full/10.1111/1752-1688.12316
	* FLUXOS-OVERLAND
		* Poster: https://www.researchgate.net/publication/333324452_Hydrodynamic_modelling_of_snowmelt_flooding_events_and_nutrient_transport_in_the_Canadian_Prairies_using_the_FLUXOS_model?channel=doi&linkId=5ce70f0a458515712ebda98b&showFulltext=true

## Branches
* master: All changes made in adesolver and adesolver_wintra have been merged into master
* developnment: branch used for development before merge with master
* adesolver: Adection-Dispersion-Reaction equation solver was converted to C++ and debugged
* adesolver_wintra: the wintra algorithm was added

## Compiling
* cmake: CMakeList is provided
* Library dependencies: Armadillo 
* Cmake minimum version: 3.10

<!-- ## Execution (and input files and folder needed) -->
## Execution
* Create a folder with name "Results" inside the working directory where the input files and fluxos are
* input files (see examle in Working_example folder)
	* master input file: e.g., modset
	* DEM file (Esri ASCII-format raster with headers removed ->  this will be fixed soon)
	* DEM of the basin (sub-set of the main DEM file for FLUXOS to know where the boundaries of the basin are)
	* Snowmelt timeseries (time,mm/day)
* to execute: ./fluxos_cpp "argument_1" (where "argument_1" is the mater input file)

<!-- ## Visualization of results (stored inside "Results" folder) -->
## Visualization
* For visualization of output stored in "Results" folder
* visit: https://wci.llnl.gov/simulation/computer-codes/visit/

<!-- ## Supporting scripts (post-processing) -->
## Supporting Scripts
* Used for post-processing
* "fluxos_python" folder

## Working Example
* "Working_example" folder


