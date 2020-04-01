
// Copyright 1992: Cornel Beffa and Roland Faeh
// Copyright 2013: Kashif Shaad and Diogo Costa
// Copyright 2019, Diogo Costa

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) an_col later version.
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

int main(int argc, char* argv[]) 
{   
    unsigned int n_rowl, n_coll, it = 0;
    unsigned int a, irow, icol, print_step, print_next, qmelt_rowi, timstart;
    double c0,v0,u0,hp, hpall, qmelti,ks_input; 
    bool outwritestatus;

    std::chrono::duration<double> elapsed_seconds;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    std::string coment_sim_str;

    std::string modset_flname (argv[1]);
    std::string dirpath = SplitFilename (modset_flname);
   
    // Create/Open log file
    std::ofstream logFLUXOSfile (dirpath + "/fluxos_run.log");
    std::cout << "FLUXOS"  << std::endl;
    logFLUXOSfile << "FLUXOS \n";
    logFLUXOSfile << "\n-----------------------------------------------\n" << std::endl;
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::cout << "Simulation started... " << std::ctime(&start_time)  << std::endl;
    logFLUXOSfile << "Simulation started... " << std::ctime(&start_time);
        
     // Get the size of the domain (nrow and ncol)
    get_domain_size(&n_rowl,&n_coll, modset_flname, dirpath, logFLUXOSfile);
     // Input the duration of the simulation
        
    // Initiate variables on the heap
    GlobVar ds(n_rowl+2,n_coll+2); 
    
    // input/read data
    ds.cfl = 1; // Courant condition
    // ds.dxy = 3; // grid size (structure grid) - it will actually come from DEM
    ds.ntim = 0;// maximum time step (seconds)
    //kapa = -2.    // /  -2=1.Ord ; -1=2.Ord   // KOMISCH, DASS REAL/INTEGER ->schauen bei Rolands Dateien
    //betas = 2. // Chezy (parameter)
    //ksfirow = 0.2 // Chezy (rougness) -> NEEDs to be converted into a vector with data for all cells
    ds.cvdef = 0.07; // for turbulent stress calc
    ds.nuem = 1.793e-6; // molecular dynamic viscosity (for turbulent stress calc)
    //print_step = 3600; // in seconds
    // timstart = 558000; // start of the simulation
        
    // read model set up
    read_modset(ds,modset_flname,dirpath,&print_step,&ks_input,logFLUXOSfile);
    
    //std::cout << "Simulation purpose (write comment):  ";
    //std::cin >> coment_sim_str;
    logFLUXOSfile << "Simulation: " + ds.sim_purp + "\n\n";
    
    // Request user input
    //std::cout << "Print step (s) = ";
    //std::cin >> print_step;
    logFLUXOSfile << "Print step (s) = " + std::to_string(print_step) + "\n";
    
    ds.n_row = ds.m_row - 2;
    ds.n_col = ds.m_col - 2;
    
    ds.D_coef = 0.01;
    
    //std::cout << "Roughness height (m) = ";
    //std::cin >> ks_input;
    logFLUXOSfile << "Roughness height (m) = " + std::to_string(ks_input) + "\n";
    
    //std::cout << "Cell size (m) = ";
    //std::cin >> ds.dxy;
    logFLUXOSfile << "Cell size (m) = " + std::to_string(ds.dxy) + "\n";
    
    ds.arbase = ds.dxy * ds.dxy;
    read_geo(ds,ks_input,logFLUXOSfile); // DEM
    
    //std::cout << "Increment to basin margins (m) = ";
    //std::cin >> zbinc;
    //logFLUXOSfile << "Increment to basin margins (m) = " + std::to_string(zbinc) + "\n";
    ds.ntim = read_load(ds,logFLUXOSfile); // snowmelt load
    
    //std::cout << "Simulation time (days) (Snowmelt input duration = " + std::to_string(ds.ntim/(3600*24)) + " days) = ";
    //std::cin >> ntim_days;
    //ds.ntim = ntim_days * 3600 * 24;
    logFLUXOSfile << "Simulation time (days) = " + std::to_string(ds.ntim) + " (= " + std::to_string(ds.ntim) + " sec)";
   
    // Input the soil nutrient release rate
    //std::cout << "Soil release rate (1/hour) = ";
    //std::cin >> ds.soil_release_rate;
    logFLUXOSfile << "\nSoil release rate (1/hour) = " + std::to_string(ds.soil_release_rate);
    
    // Input the soil background concentration
    //std::cout << "Soil initial background mass available for release to runoff (g) (0.txt points will be overwritten) = ";
    //std::cin >> ds.soil_conc_bckgrd;
    logFLUXOSfile << "\nSoil initial background mass available for release to runoff (g) (0.txt points will be overwritten) = " + std::to_string(ds.soil_conc_bckgrd);
    
    //std::cout << "SWE max (cm) = ";
    //std::cin >> ds.SWEmax;
    logFLUXOSfile << "\nSWE max (cm) = " + std::to_string(ds.SWEmax);
    ds.SWEmax = ds.SWEmax/100;
    //std::cout << "SWE std (cm) = ";
    //std::cin >> ds.SWEstd;
    logFLUXOSfile << "\nSWE std (cm) = " + std::to_string(ds.SWEstd) + "\n";
    ds.SWEstd = ds.SWEstd/100;
    
    timstart = initiation(ds,logFLUXOSfile);
    
    // INITIATION
    ds.hdry = (*ds.ks).at(1,1);  // temporary but basically saying that nothing will move until it reaches roughness height
        
    print_next = timstart;
    ds.tim = timstart;
    //write_results(ds,std::round(print_next));
    
    print_next = print_next + print_step;
        
    std::cout << "-----------------------------------------------\n" << std::endl;
    logFLUXOSfile << "\n-----------------------------------------------\n" << std::endl;
    
    // TIME LOOP
    while(ds.tim <= ds.ntim) 
    {              
        ds.dtfl=9.e10;
        hpall = 0.0f;
        
        // SPACE LOOP
        for(icol=1;icol<=ds.n_col;icol++)
        {
            for(irow=1;irow<=ds.n_row;irow++)
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
          
       
        // Qmelt load
        for (a=0;a<=(*ds.qmelt).col(0).n_elem;a++){
            qmelt_rowi = a;
            if ((*ds.qmelt).at(a,0) > ds.tim){       
                break;
            }
        }
        
        qmelti = (*ds.qmelt).at(qmelt_rowi,1)/(1000.*3600.*24.)*ds.dtfl;
        ds.qmelv_inc += qmelti;
         for(icol=1;icol<=ds.n_col;icol++)
        {
            for(irow=1;irow<=ds.n_row;irow++)
            {
                if (std::abs((*ds.zb).at(irow,icol)) != 99999)
                {
                    hp = std::max((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0); // adesolver hp before adding snowmelt  
                    (*ds.z).at(irow,icol) = (*ds.z).at(irow,icol) + qmelti;   
                    (*ds.h)(irow,icol)=std::max((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0);
                    if ((*ds.h)(irow,icol) <= ds.hdry)
                    {
                        (*ds.ldry).at(irow,icol)=0.0f;
                    }
                    (*ds.h0)(irow,icol) = (*ds.h)(irow,icol);
                    if (hp!=0.)
                    {          
                        (*ds.conc_SW)(irow,icol)=((*ds.conc_SW)(irow,icol)*hp+qmelti*0)/((*ds.h)(irow,icol)); //adesolver (adjustment for snowmelt)       
                    }
                }
            }
         }
                
        // FLOW SOLVERS
        if (hpall!=0) 
        {
            it++;
            hydrodynamics_calc(ds);
            adesolver_calc(ds, it);
            wintrasolver_calc(ds);
        }
        
        // PRINT RESULTS
        if (ds.tim>=print_next)
        {
             end = std::chrono::system_clock::now();
             elapsed_seconds = end-start;
             
              outwritestatus = write_results(ds,std::round(print_next),print_step,elapsed_seconds);
             
               
            if(outwritestatus == true) 
            {
                std::cout << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(print_step/60) << " || Time elapsed (min): " << elapsed_seconds.count()/60 << std::endl;
                logFLUXOSfile << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(print_step/60) << " || Time elapsed (min): " << std::to_string(elapsed_seconds.count()/60) + "\n";
                print_next = print_next + print_step;
                start = std::chrono::system_clock::now();

            } else
            {
                std::cout << "Problem when saving the results:" + print_next << std::endl;
                logFLUXOSfile << "Problem when saving the results:" + print_next;
                return 0;
            }
             
         }
    }
    
    // Simulation complete
    std::cout << "-----------------------------------------------" << std::endl;
    logFLUXOSfile << "\n-----------------------------------------------" << std::endl;
    std::cout << "Simulation complete (" << std::chrono::system_clock::now << ")"  << std::endl;
    logFLUXOSfile << "Simulation complete (" << std::chrono::system_clock::now;
    logFLUXOSfile.close(); 
      
    return 0;
}

