

// This program, openWQ, is free software: you can redistribute it and/or modify
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

#ifndef OPENWQ_HYDROLINK_INCLUDED
#define OPENWQ_HYDROLINK_INCLUDED


#include "openwq/src/OpenWQ_couplercalls.h"
#include "openwq/src/OpenWQ_global.h"
#include "openwq/src/OpenWQ_readjson.h"
#include "openwq/src/OpenWQ_initiate.h"
#include "openwq/src/OpenWQ_chem.h"
#include "openwq/src/OpenWQ_watertransp.h"
#include "openwq/src/OpenWQ_extwatflux_ss.h"
#include "openwq/src/OpenWQ_units.h"
#include "openwq/src/OpenWQ_solver.h"
#include "openwq/src/OpenWQ_output.h"


class ClassWQ_OpenWQ
{

    public:

    // Variables from CRHM
    const float *hru_t; // has to be converted to soil temperatures. How?
    const float *hru_area; //
    const float *SWE;
    float *SWE_conc;
    float **SWE_conc_lay;
    const float *soil_runoff;
    float *soil_runoff_cWQ;
    float **soil_runoff_cWQ_lay; 
    const float *soil_ssr;
    float *soil_ssr_conc;
    float **soil_ssr_conc_lay; 
    const float *soil_lower;
    float *soil_lower_conc;
    float **soil_lower_conc_lay;
    const float *soil_rechr;
    float *conc_soil_rechr;
    float **conc_soil_rechr_lay;
    float *surfsoil_solub_mWQ;
    float **surfsoil_solub_mWQ_lay; 
    float *conc_soil_lower;   // concentration of organic nitrogen *** from soilstate
    float **conc_soil_lower_lay;
    const float *soil_moist;
    float *Sd;
    float *Sd_conc;
    float **Sd_conc_lay;
    float *gw;
    float *gw_conc;
    float **gw_conc_lay;
    const float *soil_rechr_max;

    //ClassWQ_OpenWQ* klone(string name) const;

    void decl(
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
        OpenWQ_chem& OpenWQ_chem,                   // biochemistry modules
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
        OpenWQ_output& OpenWQ_output,                // output modules (needed for console/logfile)
        unsigned long num_hru);

    void run(
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
        OpenWQ_chem& OpenWQ_chem,                   // biochemistry modules
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
        OpenWQ_solver& OpenWQ_solver,                // solver module
        OpenWQ_output& OpenWQ_output);               // output modules


};


#endif