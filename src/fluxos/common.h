
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
#include <sys/stat.h>

#include "jnlohmann/json.hpp"
using json = nlohmann::json;

#include "../openwq/OpenWQ_hydrolink.h"


std::string SplitFilename (
    const std::string& str);

int getIntNumberFromString(
    std::string s);

double getFloatNumberFromString(
    std::string s);

// read file names in Results directory
int findLastStep(
    const char *path);

bool get_domain_size(unsigned int *rown, 
    unsigned int *coln, 
    json master_MODSET_local,
    std::ofstream& logFLUXOSfile);

bool add_meteo(
    GlobVar& ds,
    openwq_hydrolink& openwq_hydrolink,
    OpenWQ_couplercalls& OpenWQ_couplercalls,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
    OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
    OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_readjson& OpenWQ_readjson,            // read json files
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
    OpenWQ_watertransp& OpenWQ_watertransp,      // transport modules
    OpenWQ_chem& OpenWQ_chem,                    // biochemistry modules
    OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
    OpenWQ_solver& OpenWQ_solver,                // solver module
    OpenWQ_output& OpenWQ_output,
    int nchem);
    
bool add_inflow(
    GlobVar& ds,
    openwq_hydrolink& openwq_hydrolink,
    OpenWQ_couplercalls& OpenWQ_couplercalls,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
    OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
    OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_readjson& OpenWQ_readjson,            // read json files
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
    OpenWQ_watertransp& OpenWQ_watertransp,      // transport modules
    OpenWQ_chem& OpenWQ_chem,                    // biochemistry modules
    OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
    OpenWQ_solver& OpenWQ_solver,                // solver module
    OpenWQ_output& OpenWQ_output,
    int nchem);

void check_mkdir(
    std::string &dirname);

#endif