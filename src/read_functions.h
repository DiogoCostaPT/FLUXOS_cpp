
#ifndef READ_FUNCTIONSH_INCLUDED
#define READ_FUNCTIONSH_INCLUDED

#include "GlobVar.h"
#include <iostream>

void read_modset(GlobVar& ds, const std::string& filename, 
                const std::string& pathfile, 
                unsigned int *print_step, double *ks_input,
                std::ofstream& logFLUXOSfile);

void read_geo(GlobVar& ds,double ks_input,std::ofstream& logFLUXOSfile);

float read_load(GlobVar& ds,std::ofstream& logFLUXOSfile);

#endif