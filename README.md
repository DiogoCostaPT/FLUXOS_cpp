# FLUXOS-SnoWHAT
* Soure code for the FLUXOS-SnoWHAT model (FLUXOS-Snowmelt * Watershed Hydrodynamic And Transport model).
The original code (named as FLUXOS) was written in FORTRAN and consisted of the coupling of 2dmb, +QeS2, MODFLOW and MT3DMS: https://www.sciencedirect.com/science/article/pii/S1364815216306193?via%3Dihub
* Modifications:
	* Converted to C++
	* Uses: Armadillo template-based C++ library for linear algebra 
	* Removed: MODFLOW and MT3DMS (currently there is no baseflow)
	* WINTRA algorithm was integrated for calculation of runoff-soil interactions and nutrient release: https://onlinelibrary.wiley.com/doi/full/10.1002/hyp.11346

# Branches
* master: All changes made in adesolver and adesolver_wintra have been merged into master
* adesolver: Adection-Dispersion-Reaction equation solver was converted to C++ and debugged
* adesolver_wintra: the wintra algorithm was added

# Compiling
* cmake: CMakeList is provided
* Library dependencies: Armadillo 
* Cmake minimum version: 3.10

# Execution
* to execute: ./fluxos_cpp
* input files (see Working Example folder)
	* main input file: modset.fluxos
	* DEM available
	* DEM of the basin (sub-set of the main DEM that FLUXOS used to know the boundaries of the basin)
	* Snowmelt timeseries

# Visualization of results (stored inside "Results" folder)
* visit: https://wci.llnl.gov/simulation/computer-codes/visit/

# Supporting scripts (post-processing)
* "fluxos_python" folder

# Working example
* "Working_example" folder


