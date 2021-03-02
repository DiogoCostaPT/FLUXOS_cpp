
#ifndef READ_FUNCTIONSH_INCLUDED
#define READ_FUNCTIONSH_INCLUDED

#include "GlobVar.h"
#include <iostream>

bool read_modset(GlobVar& ds, const std::string& filename, 
                const std::string& pathfile, 
                unsigned int *print_step, double *ks_input,
                std::ofstream& logFLUXOSfile);

bool read_geo(GlobVar& ds,double ks_input,std::ofstream& logFLUXOSfile);

float read_meteo(GlobVar& ds,std::ofstream& logFLUXOSfile);

float read_inflow(GlobVar& ds,std::ofstream& logFLUXOSfile);

#endif