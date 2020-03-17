
#ifndef READ_FUNCTIONSH_INCLUDED
#define READ_FUNCTIONSH_INCLUDED

#include "GlobVar.h"

void read_modset(GlobVar& ds, unsigned int *print_step, double *ks_input,std::ofstream& logFLUXOSfile);

void read_geo(GlobVar& ds,double ks_input,std::ofstream& logFLUXOSfile);

float read_load(GlobVar& ds,std::ofstream& logFLUXOSfile);

#endif