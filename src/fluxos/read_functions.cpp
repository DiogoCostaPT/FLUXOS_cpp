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
    double *ks_input,
    std::ofstream& logFLUXOSfile)
{
    
    std::string str, msg;
    bool errflag = false;
    
    try{

        // Compulsory data
        ds.dem_file = ds.master_MODSET["DEM_FILE"];
        ds.output_folder = ds.master_MODSET["OUTPUT"]["OUTPUT_FOLDER"];
        ds.print_step = ds.master_MODSET["OUTPUT"]["PRINT_STEP"];
        ds.h_min_print = ds.master_MODSET["OUTPUT"]["H_MIN_TO_PRINT"];
        *ks_input = ds.master_MODSET["ROUGNESS_HEIGHT"];
        ds.soil_release_rate = ds.master_MODSET["SOIL_RELEASE_RATE"];
        ds.soil_conc_bckgrd = ds.master_MODSET["SOIL_CONC_BACKGROUND"];
        ds.SWEstd = ds.master_MODSET["SWE_STD"];
        ds.SWEmax = ds.master_MODSET["SWE_MAX"];

        // Modules
        ds.openwq = ds.master_MODSET["EXTERNAL_MODULES"]["OPENWQ"]["STATUS"];
        if (ds.openwq == true){
            ds.openwq_masterfile = ds.master_MODSET["EXTERNAL_MODULES"]["OPENWQ"]["MASTERFILE_DIR"];
        }
        // Only requires one of these forcing files to be provided
        auto exist_meteo = ds.master_MODSET.find("METEO_FILE");
        auto exist_inflow = ds.master_MODSET.find("INFLOW_FILE");

        if (exist_meteo != ds.master_MODSET.end() || 
            exist_inflow != ds.master_MODSET.end()){

            if (exist_meteo != ds.master_MODSET.end())
                ds.meteo_file = ds.master_MODSET["METEO_FILE"];

            if (exist_inflow != ds.master_MODSET.end()){
                ds.inflow_file = ds.master_MODSET["INFLOW_FILE"]["FILENAME"];

                ds.inflow_ycoord = ds.master_MODSET["INFLOW_FILE"]["DISCHARGE_LOCATION"]["Y_COORDINATE"];
                ds.inflow_xcoord = ds.master_MODSET["INFLOW_FILE"]["DISCHARGE_LOCATION"]["X_COORDINATE"];

            }
            
            msg = "Successful loading of master input file: " + filename;
        }
        else{
            msg = "METEO_FILE or INFLOW_FILE must be provided in master file: " + filename;
            errflag = true;
        }               

    }catch(json::type_error& e){ 

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
    (ds.XLLCORNER) = getFloatNumberFromString(line);
    getline (myfile,line); // YLLCORNER
    (ds.YLLCORNER) = getFloatNumberFromString(line);
    getline (myfile,line); // CELLSIZE
    (ds.dxy) =  getIntNumberFromString(line);
    getline (myfile,line); // NODATA_VALUE
    (ds.NODATA_VALUE) =  getIntNumberFromString(line);


    // Read DEM from file and store in ds.zb
    if (myfile.is_open()) //if the file is open
    {
        irow = ds.NROWS;
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
            irow--;
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

    // Determine discharge row and col
    if (ds.inflow_ycoord != NULL || ds.inflow_xcoord != NULL){
        ds.inflow_nrow = std::round(ds.inflow_ycoord - ds.YLLCORNER)/ds.dxy; 
        ds.inflow_ncol = std::round(ds.inflow_xcoord - ds.XLLCORNER)/ds.dxy;
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

            // For WINTRA 
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

    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load(ds.inflow_file,arma::csv_ascii);
    if(flstatusQ == true) {

        for(a=1;a<filedataQ.col(1).n_elem;a++){ // a == 1 because the first line is the header
            tinflows = filedataQ(a,0);  // t melt seconds
            vinflow = filedataQ(a,1);  // value of melt
            (*ds.inflow).at(a-1,0) = tinflows;  
            (*ds.inflow).at(a-1,1) = vinflow;
        
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