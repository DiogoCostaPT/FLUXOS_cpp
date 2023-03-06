
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

#include <iostream>
#include <fstream>
#include <math.h>
#include <armadillo>
#include <string>
#include <memory> 
#include <chrono>
#include <ctime>  

#include <vector> 
#include <dirent.h>
#include <sys/types.h>

#include "jnlohmann/json.hpp"
using json = nlohmann::json;

#include "common.h"
#include "GlobVar.h"
#include "read_functions.h"
#include "initiate.h"
#include "solver_drydomain.h"
#include "solver_wetdomain.h"
#include "hydrodynamics_calc.h"
#include "ADEsolver_calc.h"
#include "WINTRAsolver_calc.h"
#include "write_results.h"

// Openwq objects
#include "../openwq/OpenWQ_hydrolink.h"

OpenWQ_couplercalls OpenWQ_couplercalls;
OpenWQ_hostModelconfig OpenWQ_hostModelconfig;
OpenWQ_json OpenWQ_json;                    // create OpenWQ_json object
OpenWQ_wqconfig OpenWQ_wqconfig;            // create OpenWQ_wqconfig object
OpenWQ_units OpenWQ_units;                  // functions for unit conversion
OpenWQ_utils OpenWQ_utils;
OpenWQ_readjson OpenWQ_readjson;            // read json files
int num_HydroComp = 1; // number of compartments to link openWQ (see details in OpenWQ_hydrolink.cpp) 
int num_EWF = 2;
OpenWQ_vars OpenWQ_vars(num_HydroComp, num_EWF);
OpenWQ_initiate OpenWQ_initiate;            // initiate modules
OpenWQ_watertransp OpenWQ_watertransp;      // transport modules
OpenWQ_chem OpenWQ_chem;                    // biochemistry modules
OpenWQ_extwatflux_ss OpenWQ_extwatflux_ss;        // sink and source modules
OpenWQ_solver OpenWQ_solver;                // solver module
OpenWQ_output OpenWQ_output;                // output module
openwq_hydrolink openwq_hydrolink;

// SW id in openwq
int openwq_cmp_sw_id = 0;

int main(int argc, char* argv[]) 
{   
    unsigned int NROWSl, NCOLSl, it = 0;
    unsigned int irow, icol, print_next, timstart;
    int ntim_meteo, ntim_inflow;
    double c0,hp,v0,u0, hpall, ks_input; 
    bool outwritestatus;
    int nchem;
    std::vector<unsigned int> chem_mobile;

    bool errflag = false;

    std::chrono::duration<double> elapsed_seconds;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    std::string coment_sim_str;

    std::string modset_flname (argv[1]);
    std::string dirpath = SplitFilename (modset_flname);
    
    // #######################################################
    // create/Open Log file
    // #######################################################
    std::ofstream logFLUXOSfile (dirpath + "/fluxos_run.log");
    std::cout << "FLUXOS"  << std::endl;
    logFLUXOSfile << "FLUXOS \n";
    logFLUXOSfile << "\n-----------------------------------------------\n" << std::endl;
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::cout << "Simulation started... " << std::ctime(&start_time)  << std::endl;
    logFLUXOSfile << "Simulation started... " << std::ctime(&start_time);
    
    // #######################################################
    // Read modset MASTER FILE (JSON)
    // #######################################################
    json master_MODSET_local;
    std::ifstream i(modset_flname);
    master_MODSET_local = json::parse(i,
        /* callback */ nullptr,
        /* allow exceptions */ false,
        /* skip_comments */ true);
    
    // #######################################################
    // Get the size of the domain (nrow and ncol)
    // #######################################################
    errflag = get_domain_size(
        &NROWSl,
        &NCOLSl,
        master_MODSET_local,
        logFLUXOSfile);
    if (errflag)
        exit(EXIT_FAILURE);
    
    // #######################################################
    // Initiate variables on the heap
    // #######################################################
    GlobVar ds(
        NROWSl+2,
        NCOLSl+2); 
    
    // #######################################################
    // Input/read data
    // #######################################################
    ds.cfl = 1; // Courant condition
    ds.ntim = 0;// maximum time step (seconds)
    ds.cvdef = 0.07; // for turbulent stress calc
    ds.nuem = 1.793e-6; // molecular dynamic viscosity (for turbulent stress calc)

    // #######################################################
    // Save master_MODSET_local in GlobVar
    // #######################################################
    ds.master_MODSET = master_MODSET_local;
    
    // #######################################################
    // read model set up
    // #######################################################
    errflag = read_modset(
        ds,
        modset_flname, 
        &ks_input,
        logFLUXOSfile);
    if (errflag)
        exit(EXIT_FAILURE);

    // #######################################################
    // Provide info to console
    // #######################################################
    logFLUXOSfile << "Simulation: " + ds.sim_purp + "\n\n";
    logFLUXOSfile << "Print step (s) = " + std::to_string(ds.print_step) + "\n";
    logFLUXOSfile << "Roughness height (m) = " + std::to_string(ks_input) + "\n";
    logFLUXOSfile << "Cell size (m) = " + std::to_string(ds.dxy) + "\n";

    // #######################################################
    // Set size of domain and other info
    // #######################################################
    ds.NROWS = ds.MROWS - 2;
    ds.NCOLS = ds.MCOLS - 2;
    
    // #######################################################
    // Read DEM
    // #######################################################    
    errflag = read_geo(
        ds,
        ks_input,
        logFLUXOSfile); // DEM
    if (errflag)
        exit(EXIT_FAILURE);
    ds.arbase = ds.dxy * ds.dxy;
    
    // #######################################################
    // Read forcing: Meteo and inflow files
    // #######################################################
    ntim_meteo = read_meteo(
        ds,
        logFLUXOSfile); //  load
    ntim_inflow = read_inflow(
        ds,
        logFLUXOSfile); //  load


    // #######################################################
    // Provide simulation duraction to console
    // #######################################################
    ds.ntim = std::max(ntim_meteo,ntim_inflow); // get the max ntim

    // #######################################################
    // Provide additional information to console
    // #######################################################
    logFLUXOSfile << "Simulation time (days) = " + std::to_string(ds.ntim) + " (= " + std::to_string(ds.ntim) + " sec)";
    logFLUXOSfile << "\nSoil release rate (1/hour) = " + std::to_string(ds.soil_release_rate);
    logFLUXOSfile << "\nSoil initial background mass available for release to runoff (g) (0.txt points will be overwritten) = " + std::to_string(ds.soil_conc_bckgrd);
    logFLUXOSfile << "\nSWE max (cm) = " + std::to_string(ds.SWEmax);
    logFLUXOSfile << "\nSWE std (cm) = " + std::to_string(ds.SWEstd) + "\n";

    // #######################################################
    // Modules 
    // #######################################################
    // ade_solver
    if (ds.ade_solver == true){
        // Message that openwq has been activated
        logFLUXOSfile << "\n > ADE_solver activated";
    }else{
        // otherwise disable wintra and openwq
        ds.wintra = false;
        ds.openwq = false;
    }
    
    if (ds.openwq == true){

        // Message that openwq has been activated
        logFLUXOSfile << "\n > OpenWQ activated";

        // Call openwq_decl
        openwq_hydrolink.openwq_decl(
            OpenWQ_couplercalls,     // Class with all call from coupler
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
            OpenWQ_output,
            ds.openwq_masterfile,
            ds.MROWS,
            ds.MCOLS
        );

        nchem = OpenWQ_wqconfig.BGC_general_num_chem;
        chem_mobile = OpenWQ_wqconfig.BGC_general_mobile_species;

    }else{
        nchem = 1;
        std::vector<unsigned int> vecTemp(1);
        chem_mobile = vecTemp;
        chem_mobile[0] = 0;
    }


    // #######################################################
    // Initiate
    // #######################################################
    timstart = initiation(
        ds,
        nchem,
        logFLUXOSfile);
    ds.hdry = (*ds.ks).at(1,1);  // temporary but basically saying that nothing will move until it reaches roughness height
    print_next = timstart;  
    print_next = print_next + ds.print_step;
    if(ds.wintra==true){ds.SWEstd = ds.SWEstd/100;}

    std::cout << "-----------------------------------------------\n" << std::endl;
    logFLUXOSfile << "\n-----------------------------------------------\n" << std::endl;
    
    
    // #######################################################
    // Courant Condition: determine maximum time step for numerical stabilitity
    // #######################################################
    while(ds.tim <= ds.ntim) 
    {   
        // Reset vars
        ds.dtfl=9.e10;
        hpall = 0.0f;

        // #######################################################
        // OpenWQ run_time_start
        // #######################################################
        if (ds.openwq == true){

            openwq_hydrolink.openwq_time_start(
                OpenWQ_couplercalls,
                OpenWQ_hostModelconfig,
                OpenWQ_json,                    // create OpenWQ_json object
                OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
                OpenWQ_units,                  // functions for unit conversion
                OpenWQ_utils,
                OpenWQ_readjson,            // read json files
                OpenWQ_vars,
                OpenWQ_initiate,            // initiate modules
                OpenWQ_watertransp,      // transport modules
                OpenWQ_chem,                    // biochemistry modules
                OpenWQ_extwatflux_ss,        // sink and source modules)
                OpenWQ_solver,                // solver module
                OpenWQ_output,                // output modules
                ds);
        }
        
        for(icol=1;icol<=ds.NCOLS;icol++)
        {
            for(irow=1;irow<=ds.NROWS;irow++)
            {
                hp = (*ds.h).at(irow,icol);
                (*ds.h0)(irow,icol) = hp; // adesolver
                (*ds.ldry_prev).at(irow,icol) = (*ds.ldry).at(irow,icol); // adesolver

                if(hp>ds.hdry) 
                {
                    (*ds.ldry).at(irow,icol)=0.0f;
                    hp=std::fmax((*ds.h).at(irow,icol),ds.hdry);
                    hpall = std::fmax(hpall,(*ds.h).at(irow,icol));
                    c0=sqrt(ds.gacc*(*ds.h).at(irow,icol));
                    u0=std::fmax(.000001,fabs((*ds.qx).at(irow,icol)/hp));
                    v0=std::fmax(.000001,fabs((*ds.qy).at(irow,icol)/hp));
                    ds.dtfl=fmin(fmin(ds.cfl*ds.dxy/(u0+c0),ds.cfl*ds.dxy/(v0+c0)),ds.dtfl);
                }else 
                {
                    (*ds.ldry).at(irow,icol)=1.0f;
                }
                ds.dtfl=fmin(print_next - ds.tim, ds.dtfl);
            }
        }
                
        ds.tim = ds.tim + ds.dtfl;
          
        // #######################################################
        // Add forcing: meteo and inflow
        // #######################################################
        errflag = add_meteo(
            ds,
            openwq_hydrolink,
            OpenWQ_couplercalls,
            OpenWQ_hostModelconfig,
            OpenWQ_json,                    // create OpenWQ_json object
            OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
            OpenWQ_units,                  // functions for unit conversion
            OpenWQ_utils,
            OpenWQ_readjson,            // read json files
            OpenWQ_vars,
            OpenWQ_initiate,            // initiate modules
            OpenWQ_watertransp,      // transport modules
            OpenWQ_chem,                    // biochemistry modules
            OpenWQ_extwatflux_ss,        // sink and source modules)
            OpenWQ_solver,                // solver module
            OpenWQ_output,
            nchem);

        if (errflag) exit(EXIT_FAILURE);

        errflag = add_inflow(
            ds,
            openwq_hydrolink,
            OpenWQ_couplercalls,
            OpenWQ_hostModelconfig,
            OpenWQ_json,                    // create OpenWQ_json object
            OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
            OpenWQ_units,                  // functions for unit conversion
            OpenWQ_utils,
            OpenWQ_readjson,            // read json files
            OpenWQ_vars,
            OpenWQ_initiate,            // initiate modules
            OpenWQ_watertransp,      // transport modules
            OpenWQ_chem,                    // biochemistry modules
            OpenWQ_extwatflux_ss,        // sink and source modules)
            OpenWQ_solver,                // solver module
            OpenWQ_output,
            nchem);

        if (errflag) exit(EXIT_FAILURE);

        // Return OpenWQ conc cubes to Fluxos Mats for transport calculation
        // OpenWQ express as mass => needs to be converted to FLUXOS conc that is express
        // as conc
        if (ds.openwq == true){
            for(int ichem=0;ichem<nchem;ichem++){
                (*ds.conc_SW)[ichem] = (*OpenWQ_vars.chemass)(openwq_cmp_sw_id)(ichem).slice(0)
                        / (*ds.h * ds.arbase);
            }
        }

        // #######################################################        
        // CALL FLOW SOLVERS
        // #######################################################
        if (hpall!=0) 
        {
            it++;

            // Dynamic wave solver
            hydrodynamics_calc(ds);

            // ADE solver (that is used also by openwq)
            if(ds.ade_solver==true){
                for(int ichem=0;ichem<chem_mobile.size();ichem++){
                    adesolver_calc(ds, it, chem_mobile[ichem]);
                };  
            }

            // Wintra module
            if(ds.wintra==true){
                for(int ichem=0;ichem<chem_mobile.size();ichem++){
                    wintrasolver_calc(ds, chem_mobile[ichem]);
                };
            }

        }

        // #######################################################
        // OpenWQ run_time_end
        // #######################################################
        if (ds.openwq == true){

            // Return Fluxos Mats for transport calculation to OpenWQ conc cubes
            // OpenWQ express as mass => needs to be converted to FLUXOS conc that is express
            // as conc
            for(int ichem=0;ichem<nchem;ichem++){

                 (*OpenWQ_vars.chemass)(openwq_cmp_sw_id)(ichem).slice(0) 
                    =  (*ds.conc_SW)[ichem] % (*ds.h * ds.arbase);

            }

            openwq_hydrolink.openwq_time_end(
                OpenWQ_couplercalls,     // Class with all call from coupler
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
                ds);
        }
        
        // #######################################################
        // PRINT RESULTS
        // #######################################################
        if (ds.tim>=print_next)
        {
            end = std::chrono::system_clock::now();
            elapsed_seconds = end-start;
            
            outwritestatus = write_results(
                ds,
                std::round(print_next),
                elapsed_seconds);
             
            if(outwritestatus == true) 
            {
                std::cout << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(ds.print_step/60) << " || Time elapsed (min): " << elapsed_seconds.count()/60 << std::endl;
                logFLUXOSfile << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(ds.print_step/60) << " || Time elapsed (min): " << std::to_string(elapsed_seconds.count()/60) + "\n";
                print_next = print_next + ds.print_step;
                start = std::chrono::system_clock::now();

            } else
            {
                std::cout << "Problem when saving the results:" + print_next << std::endl;
                logFLUXOSfile << "Problem when saving the results:" + print_next;
                return 0;
            }
             
         }

    }
    
    // #######################################################
    // Simulation complete
    // #######################################################
    std::cout << "-----------------------------------------------" << std::endl;
    logFLUXOSfile << "\n-----------------------------------------------" << std::endl;
    std::cout << "Simulation complete (" << std::chrono::system_clock::now << ")"  << std::endl;
    logFLUXOSfile << "Simulation complete (" << std::chrono::system_clock::now;
    logFLUXOSfile.close(); 
      
    return 0;
}

