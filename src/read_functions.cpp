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

#include <armadillo>
#include <memory> 
#include <iostream>
#include <string>

#include "read_functions.h"
#include "common.h"

bool read_modset(
    GlobVar& ds, 
    const std::string& filename, 
    const std::string& pathfile, 
    unsigned int *print_step, 
    double *ks_input,
    std::ofstream& logFLUXOSfile)
{
    // read_modset(ds,print_step,ks_input,zbinc,ntim_days)
    
    std::string str, msg;
    bool errflag = false;
    
    std::ifstream file(filename);
    
    int i = 0;
    while (std::getline(file, str)) 
    {
        // Continue if comment
        if(str.find("#") != std::string::npos)
            continue;

        i += 1;
        //if(str.find("COMMNET") != std::string::npos){ds.sim_purp = str.substr(8);}; // comment
        if(str.find("DEM_FILE") != std::string::npos){ds.dem_file = str.substr(9);}; // DEM ESRI-ArcGIS ascii
        if(str.find("METEO_FILE") != std::string::npos){ds.meteo_file = str.substr(11);}; // Meteo file
        if(str.find("INFLOW_FILE") != std::string::npos){ds.inflow_file = str.substr(12);}; // Meteo file
        if(str.find("PRINT_STEP") != std::string::npos){(*print_step) = std::stoi(str.substr(11));}; // print time step
        if(str.find("ROUGNESS_HEIGHT") != std::string::npos){(*ks_input) = std::stof(str.substr(16));};  // average roughness height (m)
        if(str.find("SOIL_RELEASE_RATE") != std::string::npos){(ds.soil_release_rate) = std::stof(str.substr(18));}; //  WINTRA: soil nutrient release rate
        if(str.find("SOIL_CONC_BACKGROUND") != std::string::npos){(ds.soil_conc_bckgrd) = std::stof(str.substr(21));};  // WINTRA: soil background concentration
        if(str.find("SWE_STD") != std::string::npos){(ds.SWEstd) = std::stof(str.substr(8));}; // SWE standard deviation (snow depletion curves, Kevin's paper)
        if(str.find("SWE_MAX") != std::string::npos){(ds.SWEmax) = std::stof(str.substr(8));};  // SWE standard deviation (snow depletion curves, Kevin's paper)   
    }
    file.close();
    
    if(i==8){
        msg = "Successful loading of master input file: " + filename;
    } else{
        msg = "PROBLEM loading of master input file: " + filename;
        errflag = true;
    } 
     std::cout << msg  << std::endl;
     logFLUXOSfile << msg + "\n";

     return errflag;
    
}

bool read_geo(
    GlobVar& ds,
    double ks_input,
    std::ofstream& logFLUXOSfile)
{
    int icol,irow,n;  
    float zbp,zbp_corr,zbn,zbs,zbe,zbw,temp_float;
    arma::mat filedata; 
    std::string msg, temp_str;
    bool errflag = false;

    arma::mat zb_raw; 
    zb_raw = arma::zeros<arma::mat>(ds.MROWS,ds.MCOLS);

    std::string line; //this will contain the data read from the file
    std::ifstream myfile (ds.dem_file); //opening the file.
    std::stringstream str_strm;

    getline (myfile,line); // NCOLS
    getline (myfile,line); // NROWS
    getline (myfile,line); // XLLCORNER
    getline (myfile,line); // YLLCORNER
    getline (myfile,line); // CELLSIZE
    (ds.dxy) =  getIntNumberFromString(line);
    getline (myfile,line); // NODATA_VALUE
    (ds.NODATA_VALUE) =  getIntNumberFromString(line);


    // Read DEM from file and store in ds.zb
    if (myfile.is_open()) //if the file is open
    {
        irow = 1;
        while (! myfile.eof() ) //while the end of file is NOT reached
        {
            std::stringstream str_strm;
            getline (myfile,line); //get one line from the file
            str_strm << line; //convert the string s into stringstream
            
            icol = 1;
            while(!str_strm.eof()) {
                str_strm >> temp_str; //take words into temp_str one by one
                if(std::stringstream(temp_str) >> temp_float) { //try to convert string to int
                    zb_raw.at(irow,icol) = temp_float;  // For now -99999 is set as 99999 to act like a wall using abs
                    //std::cout << std::to_string(irow) + ", "+ std::to_string(icol) + ", "+ std::to_string(zb_raw.at(irow,icol)) << std::endl;
                }
                icol++;
                temp_str = ""; //clear temp string
            }
            str_strm.str(std::string());
            irow++;
        }
        myfile.close(); //closing the file
        msg = "Successful loading of DEM file: " + ds.dem_file;
    } else{
        msg = "PROBLEM loading of DEM file: " + ds.dem_file;    
        errflag = true;   
    } 
    std::cout << msg << std::endl;
    logFLUXOSfile << msg + "\n" ;


    // Neumann condition for outer and NODATA_VALUE cells
    for(icol=1;icol<=ds.NCOLS;icol++)
    {
        for(irow=1;irow<=ds.NROWS;irow++)
        {   
            zbp = zb_raw.at(irow,icol);

            if (zbp == ds.NODATA_VALUE || zbp == 0.0f){ // check if NODATA_VALUE -> if yes, then it will behave as a weir
               
                zbn = zb_raw.at(irow+1,icol);
                zbs = zb_raw.at(irow-1,icol);
                zbe = zb_raw.at(irow,icol+1);
                zbw = zb_raw.at(irow,icol-1);
                zbp_corr = 0.0f;
                n = 0;
            
                (*ds.innerNeumannBCWeir).at(irow,icol) = 1.0f;
                if (zbn != ds.NODATA_VALUE && zbn != 0.0f && irow < ds.NROWS) {
                    zbp_corr += zbn;
                    n++;
                }
                if (zbs != ds.NODATA_VALUE && zbs != 0.0f && irow > 1) {
                    zbp_corr += zbs;
                    n++;
                }
                if (zbe != ds.NODATA_VALUE && zbe != 0.0f && icol < ds.NCOLS) {
                    zbp_corr += zbe;
                    n++;
                }
                if (zbw != ds.NODATA_VALUE && zbw != 0.0f && icol > 1) {
                    zbp_corr += zbw;
                    n++;
                }
                if (zbp_corr == 0.0f) {
                    zbp_corr = zbp;
                }else
                {
                    zbp_corr /= n;
                }
                
            }else
            {
                (*ds.innerNeumannBCWeir).at(irow,icol) = 0.0f;
                zbp_corr = std::fabs(zbp);
            }
            //zb_raw.at(irow,icol) = zbp_corr;
            (*ds.zb).at(irow,icol) = zbp_corr;  // For now -99999 is set as 99999 to act like a wall using abs
            //(*ds.z).at(irow,icol) = zbp_corr;
            (*ds.ks).at(irow,icol) = ks_input; 
        }
    }

    return errflag;

}

float read_meteo(
    GlobVar& ds,
    std::ofstream& logFLUXOSfile)
{
    unsigned int a; 
    unsigned int icol,irow;
    double tmeteo,vmeteo, tmeteo_bef = 0.0f;
    double tim;
    std::string msg;
    
    // Return if no METEO_FILE provided
    if (ds.meteo_file.empty()){
        tim = 0.0;
        return tim;
    }

    // reading qmelt 
    ds.qmelvtotal  = 0;
    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load(ds.meteo_file,arma::csv_ascii);

    if(flstatusQ == true) {
        for(a=1;a<filedataQ.col(1).n_elem;a++){ // a == 1 because the first line is the header
            tmeteo = filedataQ(a,0);  // t melt seconds
            vmeteo = filedataQ(a,1);  // value of melt
            (*ds.meteo).at(a,0) = tmeteo;  
            (*ds.meteo).at(a,1) = vmeteo;
            ds.qmelvtotal += vmeteo /(1000.*3600.*24.) * (tmeteo - tmeteo_bef);  // input in mm/day
            tmeteo_bef = tmeteo;
        }
       msg = "Successful loading of METEO file: " + ds.meteo_file;
    } else{
       msg = "PROBLEM loading of METEO file: " + ds.meteo_file;       
    } 
    std::cout << msg  << std::endl;
    logFLUXOSfile << msg + "\n";
    
    tim = tmeteo;
    return tim;
}


float read_inflow(
    GlobVar& ds,
    std::ofstream& logFLUXOSfile)
{
    unsigned int a; 
    unsigned int icol,irow;
    double tinflows,vinflow, tinflows_bef = 0.0f;
    double tim;
    std::string msg;

    // Return if no METEO_FILE provided
    if (ds.inflow_file.empty()){
        tim = 0.0;
        return tim;
    }
    
    // reading inflow
    ds.qmelvtotal  = 0;
    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load(ds.inflow_file,arma::csv_ascii);
    if(flstatusQ == true) {

        ds.ix_inflow = filedataQ(1,0);
        ds.iy_inflow = filedataQ(1,1);

        for(a=2;a<filedataQ.col(1).n_elem;a++){ // a == 1 because the first line is the header
            tinflows = filedataQ(a,0);  // t melt seconds
            vinflow = filedataQ(a,1);  // value of melt
            (*ds.inflow).at(a-2,0) = tinflows;  
            (*ds.inflow).at(a-2,1) = vinflow;
            ds.qmelvtotal += vinflow /(1000.*3600.*24.) * (tinflows - tinflows_bef);  // input in mm/day
            tinflows_bef = tinflows;
        }
       msg = "Successful loading of INFLOW file: " + ds.inflow_file;
    } else{
        msg = "PROBLEM loading of INFLOW file: " + ds.inflow_file;       
    } 
    std::cout << msg  << std::endl;
    logFLUXOSfile << msg + "\n";
    
    tim = tinflows;
    return tim;
}