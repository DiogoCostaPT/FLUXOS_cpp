// Copyright 2020, Diogo Costa, diogo.pinhodacosta@canada.ca
// This file is part of OpenWQ model.

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


#include "OpenWQ_hydrolink.h"


void openwq_hydrolink::openwq_decl(
    OpenWQ_couplercalls& OpenWQ_couplercalls,     // Class with all call from coupler
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
    OpenWQ_output& OpenWQ_output,
    std::string openwq_masterfile,
    unsigned long MROWS,
    unsigned long MCOLS)
{
    
    // Location of master file
    OpenWQ_wqconfig.OpenWQ_masterjson = openwq_masterfile;
    // Initiate output VTU file name string 
    //STASKED// std::string vtufilename;

    // Check if the call is from ::decl (allow to set OpenWQ_hostModelconfig)
    // or ::initbase (do not allow because OpenWQ_hostModelconfig has already been
    // defined)
    if (OpenWQ_hostModelconfig.HydroComp.size()==0){

        // #######################################
        // Characterize the Host model domain
        // Host model specific
        // #######################################
        OpenWQ_hostModelconfig.HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(0,"SURFACE_WATER",MROWS,MCOLS,1));

        // (add other compartments as needed)...

        // External fluxes
        // Make sure to use capital letters for external fluxes
        OpenWQ_hostModelconfig.HydroExtFlux.push_back(OpenWQ_hostModelconfig::hydroTuple(0,"INFLOW", MROWS,MCOLS,1)); // input at just one location
        OpenWQ_hostModelconfig.HydroExtFlux.push_back(OpenWQ_hostModelconfig::hydroTuple(1,"METEO", MROWS,MCOLS,1)); // input at just one location

        // Dependencies
        // to expand BGC modelling options
        //OpenWQ_hostModelconfig.HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(0,"SM", nhru,1,1));
        //OpenWQ_hostModelconfig.HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(1,"Tair_K", nhru,1,1));
        //OpenWQ_hostModelconfig.HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(2,"Tsoil_K", nhru,1,1));

        // Create Object: OpenWQ_json (Input JSON files) and wqconfig
        //OpenWQ_json OpenWQ_json;            // create OpenWQ_json object
        //OpenWQ_wqconfig OpenWQ_wqconfig;    // create OpenWQ_wqconfig object
        //OpenWQ_units OpenWQ_units;          // functions for unit conversion
        
        OpenWQ_couplercalls.InitialConfig(
            OpenWQ_hostModelconfig,
            OpenWQ_json,                    // create OpenWQ_json object
            OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
            OpenWQ_units,                  // functions for unit conversion
            OpenWQ_utils,
            OpenWQ_readjson,            // read json files
            OpenWQ_vars,
            OpenWQ_initiate,            // initiate modules
            OpenWQ_watertransp,      // transport modules
            OpenWQ_chem,                   // biochemistry modules
            OpenWQ_extwatflux_ss,        // sink and source modules)
            OpenWQ_output);




    }

    
}


void openwq_hydrolink::openwq_time_start(
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
    OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,  // sink and source modules)
    OpenWQ_solver& OpenWQ_solver,                // solver module
    OpenWQ_output& OpenWQ_output,                // output modules
    GlobVar& GlobVar_fluxos)                    
{

    // Local Variables
    unsigned int Sub_mob;   // interactive species index (for mobile species)
    unsigned long fluxos_nrows = GlobVar_fluxos.MROWS;
    unsigned long fluxos_ncols = GlobVar_fluxos.MCOLS;

    // Retrieve simulation timestamp
    // convert to OpenWQ time convention: seconds since 00:00 hours, Jan 1, 1900 UTC
    // this allows the use of a number of funtions of C++
    time_t simtime = getSimTime(
        OpenWQ_wqconfig,
        OpenWQ_units,
        GlobVar_fluxos.sim_start_time, 
        GlobVar_fluxos.tim);
    //int hru = 10;
    // int nhru = 100;
    // int hh = 1;

    /*
    // Update Dependencies to kinetic formulas (needs loop to get hydro model variables
    // that are dependencies to OpenWQ)
    //#pragma omp parallel for num_threads(OpenWQ_wqconfig.num_threads_requested)
    for (unsigned int hru=0;hru<nhru;hru++){
        (*OpenWQ_hostModelconfig.dependVar)[0](hru,0,0) = soil_rechr[hru]/soil_rechr_max[hru];  // loop needed - Save all SM data from hostmodel at time t
        (*OpenWQ_hostModelconfig.dependVar)[1](hru,0,0) = hru_t[hru];  // loop needed - Save all Taair data from hostmodel at time t      
        (*OpenWQ_hostModelconfig.dependVar)[2](hru,0,0) = hru_t[hru];   // keeping the same as Tair for now
        //(*OpenWQ_hostModelconfig.SM)(hru,0,0) = soil_rechr[hru]/soil_rechr_max[hru];  // loop needed - Save all SM data from hostmodel at time t
        //(*OpenWQ_hostModelconfig.Tair)(hru,0,0) = hru_t[hru];  // loop needed - Save all Taair data from hostmodel at time t      
        //(*OpenWQ_hostModelconfig.Tsoil)(hru,0,0) = hru_t[hru];   // keeping the same as Tair for now
    }
    */

    // Get Fluxes from Hydrological model
    // Save them to (*OpenWQ_hostModelconfig.waterVol_hydromodel)
    // Convert units to m3

    // CRHM water units are in mm
    // CRHM area units are in km2
    // So conversion of mm (water) to m3 is calculated as follows:
    // water_m3 = (water_mm / 1000) * (area_km2 * 1000000) <=>
    // water_m3 = water_mm * area_km2 * 1000
    #pragma omp parallel for private (fluxos_nrows, fluxos_ncols) num_threads(OpenWQ_wqconfig.num_threads_requested)
    for (int ir=0;ir<fluxos_nrows;ir++){
        for (int ic=0;ic<fluxos_ncols;ic++){

            // Surface water
            (*OpenWQ_hostModelconfig.waterVol_hydromodel)[0](ir,ic,0) 
                = std::fmax((*GlobVar_fluxos.h).at(ic,ir) * GlobVar_fluxos.arbase 
                    , 0.0f);

        }
       
    }

    // ##############################################
    // Calls all functions required
    // INSIDE the TIME loop
    // But BEFORE the SPACE loop is initiated
    // ##############################################
    OpenWQ_couplercalls.RunTimeLoopStart(
            OpenWQ_hostModelconfig,
            OpenWQ_json,                    // create OpenWQ_json object
            OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
            OpenWQ_units,                  // functions for unit conversion
            OpenWQ_utils,
            OpenWQ_readjson,            // read json files
            OpenWQ_vars,
            OpenWQ_initiate,            // initiate modules
            OpenWQ_watertransp,      // transport modules
            OpenWQ_chem,                   // biochemistry modules
            OpenWQ_extwatflux_ss,        // sink and source modules)
            OpenWQ_solver,
            OpenWQ_output,
            simtime);

}

void openwq_hydrolink::run_space_in(
    GlobVar& GlobVar,
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
    OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,  // sink and source modules)
    OpenWQ_solver& OpenWQ_solver,                // solver module
    OpenWQ_output& OpenWQ_output,                // output modules
    std::string source_EWF_name,
    int ix_r, int iy_r,
    double wflux_s2r){

    // Local variables
    int recipient = 0;
    int iz_r = 0;

// Retrieve simulation timestamp
    // convert to OpenWQ time convention: seconds since 00:00 hours, Jan 1, 1900 UTC
    // this allows the use of a number of funtions of C++
    time_t simtime = getSimTime(
        OpenWQ_wqconfig,
        OpenWQ_units,
        GlobVar.sim_start_time, 
        GlobVar.tim);

    // call RunSpaceStep_IN
    OpenWQ_couplercalls.RunSpaceStep_IN(
        OpenWQ_hostModelconfig,
        OpenWQ_json,
        OpenWQ_wqconfig,
        OpenWQ_units,
        OpenWQ_utils,
        OpenWQ_readjson,
        OpenWQ_vars,
        OpenWQ_initiate,
        OpenWQ_watertransp,
        OpenWQ_chem,
        OpenWQ_extwatflux_ss,
        OpenWQ_solver,
        OpenWQ_output,
        simtime,
        source_EWF_name,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r);

}

    /*
    
    // ##############################################
    // Calls all functions required
    // INSIDE the SPACE loop
    // TRANSPORT routines
    // ##############################################

    // Initiation steps that only need to be done once (at the start of the simulation)
    if (OpenWQ_hostModelconfig.interaction_step > 1){

        // ################################
        // Retried current state of state variables form CRHM
        // Convert from concentration (default unit in CRHM) to mass (default unit of OpenWQ)
        // ################################

        // iteractive compartment index
        unsigned int icmp;

        #pragma omp parallel for private(Sub_mob) num_threads(OpenWQ_wqconfig.num_threads_requested)
        for(unsigned int hru = 0; hru < nhru; ++hru) {

            for(unsigned int Sub_i = 0; Sub_i < OpenWQ_wqconfig.BGC_general_mobile_species.size(); ++Sub_i) {

                // ################################
                // Adv-Disp directly from CRHM's modules

                // Get mobile species indexes
                // Only update mass/concentrations using CRHM's water for mobile species
                Sub_mob = OpenWQ_wqconfig.BGC_general_mobile_species[Sub_i];
    
                // SWE
                (*OpenWQ_vars.d_chemass_dt_transp)(0)(Sub_mob)(hru,0,0) = 
                    (std::fmax(SWE_conc_lay[Sub_mob][hru],0.0f)
                        * (*OpenWQ_hostModelconfig.waterVol_hydromodel)[0](hru,0,0) 
                    ) - (*OpenWQ_vars.chemass)(0)(Sub_mob)(hru,0,0);

                // RUNOFF
                (*OpenWQ_vars.d_chemass_dt_transp)(1)(Sub_mob)(hru,0,0) = 
                    (std::fmax(soil_runoff_cWQ_lay[Sub_mob][hru],0.0f)
                    * (*OpenWQ_hostModelconfig.waterVol_hydromodel)[1](hru,0,0) 
                    ) - (*OpenWQ_vars.chemass)(1)(Sub_mob)(hru,0,0);

                // SSR
                (*OpenWQ_vars.d_chemass_dt_transp)(2)(Sub_mob)(hru,0,0) = 
                    (std::fmax(soil_ssr_conc_lay[Sub_mob][hru],0.0f)
                    * (*OpenWQ_hostModelconfig.waterVol_hydromodel)[2](hru,0,0) 
                    ) - (*OpenWQ_vars.chemass)(2)(Sub_mob)(hru,0,0);

                // SD
                (*OpenWQ_vars.d_chemass_dt_transp)(3)(Sub_mob)(hru,0,0) = 
                    (std::fmax(Sd_conc_lay[Sub_mob][hru],0.0f)
                    * (*OpenWQ_hostModelconfig.waterVol_hydromodel)[3](hru,0,0) 
                    ) - (*OpenWQ_vars.chemass)(3)(Sub_mob)(hru,0,0);

                
                // SOIL_RECHR
                (*OpenWQ_vars.d_chemass_dt_transp)(4)(Sub_mob)(hru,0,0) = 
                    (std::fmax(conc_soil_rechr_lay[Sub_mob][hru],0.0f)
                    * (*OpenWQ_hostModelconfig.waterVol_hydromodel)[4](hru,0,0) 
                    ) - (*OpenWQ_vars.chemass)(4)(Sub_mob)(hru,0,0);

                // SOIL LOWER
                (*OpenWQ_vars.d_chemass_dt_transp)(5)(Sub_mob)(hru,0,0) = 
                    (std::fmax(conc_soil_lower_lay[Sub_mob][hru],0.0f)
                    * (*OpenWQ_hostModelconfig.waterVol_hydromodel)[5](hru,0,0) 
                    ) - (*OpenWQ_vars.chemass)(5)(Sub_mob)(hru,0,0);

                // SURFSOIL 
                // (does not multiply by volume of water because CRHM uses mass 
                //for this specific compartment)
                (*OpenWQ_vars.d_chemass_dt_transp)(6)(Sub_mob)(hru,0,0) = 
                    (std::fmax(surfsoil_solub_mWQ_lay[Sub_mob][hru],0.0f)
                    ) - (*OpenWQ_vars.chemass)(6)(Sub_mob)(hru,0,0);
                                            
                
                // GW
                (*OpenWQ_vars.d_chemass_dt_transp)(7)(Sub_mob)(hru,0,0) = 
                    (std::fmax(gw_conc_lay[Sub_mob][hru],0.0f)
                    * (*OpenWQ_hostModelconfig.waterVol_hydromodel)[7](hru,0,0) 
                    ) - (*OpenWQ_vars.chemass)(7)(Sub_mob)(hru,0,0);
                
            }

            // ################################
            // EROSION: Called to calculation erosion 
            // Transport is calculated by native CRHM routines
            // only need to consider 

            // Runoff
            icmp = 1;
            OpenWQ_couplercalls.RunSpaceStep(
                OpenWQ_hostModelconfig,
                OpenWQ_json,                // create OpenWQ_json object
                OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
                OpenWQ_units,               // functions for unit conversion
                OpenWQ_utils,
                OpenWQ_readjson,            // read json files
                OpenWQ_vars,
                OpenWQ_initiate,            // initiate modules
                OpenWQ_watertransp,         // transport modules
                OpenWQ_chem,                // biochemistry modules
                OpenWQ_extwatflux_ss,          // sink and source modules)
                OpenWQ_solver,
                OpenWQ_output,
                simtime,                    // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
                icmp,                       // source: runoff
                hru,                        // x
                0,                          // y
                0,                          // z
                icmp,                       // receipient: runoff (the same as source because for CRHM we'll only deal with erosion, and for simplicity we'll just put it in the same compartment)
                hru,                        // x
                0,                          // y
                0,                          // z
                soil_runoff[hru] * 1000 * hru_area[hru],                        // water flux from source compartment to receipient compartment
                (*OpenWQ_hostModelconfig.waterVol_hydromodel)[icmp](hru,0,0));  // water mass of the receipient compartment
            
            // SSR
            icmp = 2;
            OpenWQ_couplercalls.RunSpaceStep(
                OpenWQ_hostModelconfig,
                OpenWQ_json,                    // create OpenWQ_json object
                OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
                OpenWQ_units,                  // functions for unit conversion
                OpenWQ_utils,
                OpenWQ_readjson,            // read json files
                OpenWQ_vars,
                OpenWQ_initiate,            // initiate modules
                OpenWQ_watertransp,      // transport modules
                OpenWQ_chem,                   // biochemistry modules
                OpenWQ_extwatflux_ss,        // sink and source modules)
                OpenWQ_solver,
                OpenWQ_output,
                simtime,                            // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1900 UTC
                icmp,      // source: runoff
                hru,    // x
                0,      // y
                0,      // z
                icmp,      // receipient: runoff (the same as source because for CRHM we'll only deal with erosion, and for simplicity we'll just put it in the same compartment)
                hru,    // x
                0,      // y
                0,      // z
                soil_ssr[hru] * 1000 * hru_area[hru],          // water flux from source compartment to receipient compartment
                (*OpenWQ_hostModelconfig.waterVol_hydromodel)[icmp](hru,0,0));   // water mass of the receipient compartment

        }

    }

    */


// Run time end
void openwq_hydrolink::openwq_time_end(
    OpenWQ_couplercalls& OpenWQ_couplercalls,     // Class with all call from coupler
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
    OpenWQ_solver& OpenWQ_solver,
    OpenWQ_output& OpenWQ_output,
    GlobVar& GlobVar_fluxos)
{

    // Retrieve simulation timestamp
    // convert to OpenWQ time convention: seconds since 00:00 hours, Jan 1, 1900 UTC
    // this allows the use of a number of funtions of C++
    time_t simtime = getSimTime(
        OpenWQ_wqconfig,
        OpenWQ_units,
        GlobVar_fluxos.sim_start_time, 
        GlobVar_fluxos.tim);

   // ##############################################
    // Calls all functions required
    // INSIDE the TIME loop
    // But AFTER the SPACE loop is initiated
    // ##############################################
   OpenWQ_couplercalls.RunTimeLoopEnd(
       OpenWQ_hostModelconfig,
            OpenWQ_json,                    // create OpenWQ_json object
            OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
            OpenWQ_units,                  // functions for unit conversion
            OpenWQ_utils,
            OpenWQ_readjson,            // read json files
            OpenWQ_vars,
            OpenWQ_initiate,            // initiate modules
            OpenWQ_watertransp,      // transport modules
            OpenWQ_chem,                   // biochemistry modules
            OpenWQ_extwatflux_ss,        // sink and source modules)
            OpenWQ_solver,
            OpenWQ_output,
            simtime);

    //}

    // ################################
    // Save the arma::variable results back to CRHM's variables (for each compartment)
    // Convert from mass (default unit in OpenWQ) to concentrations (default unit of CRHM-WQ)
    // ################################ 
    
    // Water volume conversions
        // CRHM water units are in mm
        // CRHM area units are in km2
        // So conversion of mm (water) to m3 is calculated as follows:
        // water_m3 = (water_mm / 1000) * (area_km2 * 1000000) <=>
        // water_m3 = water_mm * area_km2 * 1000

    /*
    #pragma omp parallel for private(Sub_mob) collapse(2) num_threads(OpenWQ_wqconfig.num_threads_requested)
    for(unsigned int hru = 0; hru < nhru; ++hru) {

        for(unsigned int Sub_i = 0; Sub_i < OpenWQ_wqconfig.BGC_general_mobile_species.size(); ++Sub_i) {

                // Get mobile species indexes
                // Only update mass/concentrations using CRHM's water for mobile species
                Sub_mob = OpenWQ_wqconfig.BGC_general_mobile_species[Sub_i];

            // SWE
            if ((*OpenWQ_hostModelconfig.waterVol_hydromodel)[0](hru,0,0) 
                    > OpenWQ_hostModelconfig.watervol_minlim){

                SWE_conc_lay[Sub_mob][hru] = 
                    (*OpenWQ_vars.chemass)(0)(Sub_mob)(hru,0,0)
                    / (*OpenWQ_hostModelconfig.waterVol_hydromodel)[0](hru,0,0);
            
            }else{
                SWE_conc_lay[Sub_mob][hru] = 0.0f;
            }

            // RUNOFF
            if ((*OpenWQ_hostModelconfig.waterVol_hydromodel)[1](hru,0,0) 
                    > OpenWQ_hostModelconfig.watervol_minlim){

                soil_runoff_cWQ_lay[Sub_mob][hru] = 
                    (*OpenWQ_vars.chemass)(1)(Sub_mob)(hru,0,0)
                    / (*OpenWQ_hostModelconfig.waterVol_hydromodel)[1](hru,0,0);
            
            }else{
                soil_runoff_cWQ_lay[Sub_mob][hru] = 0.0f;
            }

            // SSR
            if ((*OpenWQ_hostModelconfig.waterVol_hydromodel)[2](hru,0,0) 
                    > OpenWQ_hostModelconfig.watervol_minlim){

                soil_ssr_conc_lay[Sub_mob][hru]= 
                    (*OpenWQ_vars.chemass)(2)(Sub_mob)(hru,0,0)
                    / (*OpenWQ_hostModelconfig.waterVol_hydromodel)[2](hru,0,0);
            
            }else{
                soil_ssr_conc_lay[Sub_mob][hru] = 0.0f;
            }

            // SD
            if ((*OpenWQ_hostModelconfig.waterVol_hydromodel)[3](hru,0,0) 
                    > OpenWQ_hostModelconfig.watervol_minlim){

                Sd_conc_lay[Sub_mob][hru] = 
                    (*OpenWQ_vars.chemass)(3)(Sub_mob)(hru,0,0)
                    / (*OpenWQ_hostModelconfig.waterVol_hydromodel)[3](hru,0,0);
            
            }else{
                Sd_conc_lay[Sub_mob][hru] = 0.0f;
            }

            // SOIL_RECHR
            if ((*OpenWQ_hostModelconfig.waterVol_hydromodel)[4](hru,0,0) 
                    > OpenWQ_hostModelconfig.watervol_minlim){

                conc_soil_rechr_lay[Sub_mob][hru] = 
                    (*OpenWQ_vars.chemass)(4)(Sub_mob)(hru,0,0)
                    / (*OpenWQ_hostModelconfig.waterVol_hydromodel)[4](hru,0,0);
            
            }else{
                conc_soil_rechr_lay[Sub_mob][hru] = 0.0f;
            }

            // SOIL LOWER
            if ((*OpenWQ_hostModelconfig.waterVol_hydromodel)[5](hru,0,0) 
                    > OpenWQ_hostModelconfig.watervol_minlim){

                conc_soil_lower_lay[Sub_mob][hru] = 
                    (*OpenWQ_vars.chemass)(5)(Sub_mob)(hru,0,0)
                    / (*OpenWQ_hostModelconfig.waterVol_hydromodel)[5](hru,0,0);;
            
            }else{
                conc_soil_lower_lay[Sub_mob][hru] = 0.0f;
            }
            
            // SURFSOIL 
            // (does not divide by volume of water because CRHM uses mass 
            //for this specific compartment)
            surfsoil_solub_mWQ_lay[Sub_mob][hru] = (*OpenWQ_vars.chemass)(6)(Sub_mob)(hru,0,0);
            
            // GW
            if ((*OpenWQ_hostModelconfig.waterVol_hydromodel)[7](hru,0,0) 
                    > OpenWQ_hostModelconfig.watervol_minlim){

                gw_conc_lay[Sub_mob][hru] = 
                    (*OpenWQ_vars.chemass)(7)(Sub_mob)(hru,0,0)
                    / (*OpenWQ_hostModelconfig.waterVol_hydromodel)[7](hru,0,0);
            
            }else{
                gw_conc_lay[Sub_mob][hru] = 0.0f;
            }

        }
    }

    */

}

// Convert time str+int to time_t
time_t openwq_hydrolink::getSimTime(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_units& OpenWQ_units,
    std::string fluxos_sim_start_time_str, 
    double fluxos_time_secs){
    
    // Local variables
    time_t fluxos_time_as_time_t;

    // Convert time string to time_t
    fluxos_time_as_time_t = OpenWQ_units.convertTime_str2time_t(
            OpenWQ_wqconfig,
            fluxos_sim_start_time_str);

    // Adding the current sim seconds
    fluxos_time_as_time_t += fluxos_time_secs;

    return fluxos_time_as_time_t;

}