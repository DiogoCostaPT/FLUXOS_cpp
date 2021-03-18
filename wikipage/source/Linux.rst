Linux
==================================

1.Armadillo
	You could download the Armadillo with the apt package manager, but the version provided is really old, which is lower than our required version 9.9. So Download the package from the website 
	
	arma.sourceforge.net/download.html
	
	Extract the package and go to the armadillo folder, then run:
	
	``cmake .``
	
	``make``
		
	``sudo make install``
	
	
2.OpenMP
	You could directly install OpenMP with the command below, which meets our version requirement 
	
	``sudo apt-get install libomp-dev``
	
3.cmake
	You could download the cmake package from the website
	
	https://cmake.org/files/
	
	Please make sure the cmake version is hogher than 3.10
	
4.Compile
	There is a CmakeLists.txt in the github repo. You could directly compile the application with the command
	
	``cmake .``
	
	``make``
	
	The executable file will be generated to the Results folder.
