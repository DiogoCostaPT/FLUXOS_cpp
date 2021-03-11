
// Copyright 1992: Cornel Beffa and Roland Faeh
// Copyright 2013: Kashif Shaad and Diogo Costa
// Copyright 2019, Diogo Costa

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef COMMONH_INCLUDED
#define COMMONH_INCLUDED

#include "GlobVar.h"
#include <string.h>

std::string SplitFilename (
    const std::string& str);

int getIntNumberFromString(
    std::string s);

// read file names in Results directory
int findLastStep(
    const char *path);

bool get_domain_size(
    unsigned int *rown, 
    unsigned int *coln, 
    const std::string& filename, 
    const std::string& pathfile,
    std::ofstream& logFLUXOSfile);

bool add_meteo(
    GlobVar& ds);
bool add_inflow(
    GlobVar& ds);

#endif