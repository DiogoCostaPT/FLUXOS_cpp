
#ifndef COMMONH_INCLUDED
#define COMMONH_INCLUDED

#include "GlobVar.h"
#include <string.h>

std::string SplitFilename (const std::string& str);

int getIntNumberFromString(std::string s);

// read file names in Results directory
int findLastStep(const char *path);
bool get_domain_size(unsigned int *rown, unsigned int *coln, 
                    const std::string& filename, 
                    const std::string& pathfile, std::ofstream& logFLUXOSfile);

bool add_qmelt(GlobVar& ds);

#endif