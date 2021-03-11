
// Copyright 1992: Cornel Beffa and Roland Faeh
// Copyright 2013: Kashif Shaad and Diogo Costa
// Copyright 2019: Diogo Costa

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

#ifndef READ_FUNCTIONSH_INCLUDED
#define READ_FUNCTIONSH_INCLUDED

#include "GlobVar.h"
#include <iostream>

bool read_modset(
    GlobVar& ds, 
    const std::string& filename, 
    unsigned int *print_step, 
    double *ks_input,
    std::ofstream& logFLUXOSfile);

bool read_geo(
    GlobVar& ds,
    double ks_input,
    std::ofstream& logFLUXOSfile);

float read_meteo(
    GlobVar& ds,
    std::ofstream& logFLUXOSfile);

float read_inflow(
    GlobVar& ds,
    std::ofstream& logFLUXOSfile);

#endif