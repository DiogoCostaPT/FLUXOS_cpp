Windows
==================================

1.Armadillo
	Download the latest stable package from the website 
	
	arma.sourceforge.net/download.html
	
	After extract, copy the include folder to a convenient directory. We will need this directroy for the cmamke file.
	
.. image:: windows1.PNG

2.OpenMP
	We could build the C++ compiler with MinGW. Please install the Windows 64bit MinGW, otherwise you are likely meet problems when working with Armadillo library. Here is the link
	
	https://sourceforge.net/projects/mingw-w64/
	
	Running the setup file and set the architecture to "x86_64", also remember the installation directory.	

	Then go to the installation folder, add the ``<mingw installation path>\i686-8.1.0-posix-dwarf-rt_v6-rev0/mingw32/include`` to the windows environment path. By now, MinGW should have been configurated successfully. 
	
	MinGW has a package manager which would be convenient to install the libraries and packages. We could setup OpenMP with it.
	
	https://sourceforge.net/projects/mingw/
	
	The installation process is quite straight forward. After installation, you could be able to use "MinGW Installation Manager". In order to ise OpenMP, we need to install "mingw-pthread" shown as the figure below
	
.. image:: windows3.PNG
	
.. image:: windows2.PNG

3. cmake and compile
	Before runing the cmake file, please write the armadillo include directory to the "target_include_directorties" in the cmake file. After that, open the cmd and go to the directroy of the "CmakeLists.txt". Then
	
	``cmake -G "Unix Makefiles"``
	
	``make``
	
	The executable file will be generated to the "Results" folder
