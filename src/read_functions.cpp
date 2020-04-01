#include <armadillo>
#include <memory> 
#include <iostream>
#include <string>
#include <iostream>
#include <fstream>

#include "read_functions.h"
#include "common.h"


void read_modset(GlobVar& ds, const std::string& filename, 
                const std::string& pathfile, unsigned int *print_step, 
                double *ks_input,
                std::ofstream& logFLUXOSfile)
{
    // read_modset(ds,print_step,ks_input,zbinc,ntim_days)
    
    std::string str, msg;
    
    std::ifstream file(filename);
    
    int i = 0;
    while (std::getline(file, str)) 
    {
        i += 1;
        if(i==1){ds.sim_purp = str;}; // comment
        if(i==2){ds.dem_file = str;}; // DEM ESRI-ArcGIS ascii
        if(i==3){ds.qmelt_file = str;}; // snowmelt file
        if(i==4){(*print_step) = std::stoi(str);}; // print time step
        if(i==5){(*ks_input) = std::stof(str);};  // average roughness height (m)
        if(i==6){(ds.soil_release_rate) = std::stof(str);}; //  WINTRA: soil nutrient release rate
        if(i==7){(ds.soil_conc_bckgrd) = std::stof(str);};  // WINTRA: soil background concentration
        if(i==8){(ds.SWEstd) = std::stof(str);}; // SWE standard deviation (snow depletion curves, Kevin's paper)
        if(i==9){(ds.SWEmax) = std::stof(str);};  // SWE standard deviation (snow depletion curves, Kevin's paper)
        
    }
    file.close();
    
    if(i==9){
        msg = "Successful loading of master input file: " + filename;
    } else{
        msg = "PROBLEM loading of master input file: " + filename;

    } 
     std::cout << msg  << std::endl;
     logFLUXOSfile << msg + "\n";
    
}

void read_geo(GlobVar& ds,double ks_input,std::ofstream& logFLUXOSfile)
{
    int icol,irow;  
    double zbp,zbp_corr,zbn,zbs,zbe,zbw, temp_float;
    arma::mat filedata; 
    std::string msg, temp_str;

    arma::mat zb_raw;
    zb_raw = (*ds.zb);
    
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
        unsigned int irow = 0;
        unsigned int icol = 0;
        while (! myfile.eof() ) //while the end of file is NOT reached
        {
            getline (myfile,line); //get one line from the file
            str_strm << line; //convert the string s into stringstream
            
            while(!str_strm.eof()) {
                str_strm >> temp_str; //take words into temp_str one by one
                if(std::stringstream(temp_str) >> temp_float) { //try to convert string to int
                        zb_raw.at(irow,icol) = std::abs(temp_float);  // For now -99999 is set as 99999 to act like a wall using abs
                }
                icol++;
                temp_str = ""; //clear temp string
            }
            irow++;
        }
        myfile.close(); //closing the file
        msg = "Successful loading of DEM file: " + ds.dem_file;
    } else{
        msg = "PROBLEM loading of DEM file: " + ds.dem_file;       
    } 
    std::cout << msg << std::endl;
    logFLUXOSfile << msg + "\n" ;


    // Neumann condition for outer and NODATA_VALUE cells
    for(icol=1;icol<=ds.NCOLS;icol++)
    {
        for(irow=1;irow<=ds.NROWS;irow++)
        {   
            zbp = zb_raw.at(irow,icol);
            zbn = zb_raw.at(irow+1,icol);
            zbs = zb_raw.at(irow-1,icol);
            zbe = zb_raw.at(irow,icol+1);
            zbw = zb_raw.at(irow,icol-1);
            zbp_corr = 0.0f;
            if (zbp == ds.NODATA_VALUE){ // check if NODATA_VALUE -> if yes, then it will behave as a weir
                (*ds.innerNeumannBCWeir).at(irow,icol) = 1.0f;
                if (zbn != ds.NODATA_VALUE) zbp_corr = zbn;
                if (zbs != ds.NODATA_VALUE) zbp_corr = (zbp_corr + zbs)/2;
                if (zbe != ds.NODATA_VALUE) zbp_corr = 2/3*zbp_corr + 1/3*zbe;
                if (zbw != ds.NODATA_VALUE) zbp_corr = 3/4*zbp_corr + 1/4*zbw;
                if (zbp_corr == 0.0f) zbp_corr = zbp;
            }else
            {
                (*ds.innerNeumannBCWeir).at(irow,icol) = 0.0f;
                zbp_corr = zbp;
            }
            (*ds.zb).at(irow,icol) = std::abs(zbp_corr);  // For now -99999 is set as 99999 to act like a wall using abs
            (*ds.ks).at(irow,icol) = ks_input; 
        }
    }
}

float read_load(GlobVar& ds,std::ofstream& logFLUXOSfile)
{
    unsigned int a; 
    unsigned int icol,irow;
    double tmelts,vmelt, tmelts_bef = 0.0f;
    std::string msg;
        
    // reading qmelt 
    ds.qmelvtotal  = 0;
    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load(ds.qmelt_file,arma::csv_ascii);
    if(flstatusQ == true) {
        for(a=0;a<filedataQ.col(1).n_elem;a++){
            tmelts = filedataQ(a,0);  // t melt seconds
            vmelt = filedataQ(a,1);  // value of melt
            (*ds.qmelt).at(a,0) = tmelts;  
            (*ds.qmelt).at(a,1) = vmelt;
            ds.qmelvtotal += vmelt /(1000.*3600.*24.) * (tmelts - tmelts_bef); 
            tmelts_bef = tmelts;
        }
       msg = "Successful loading of Qmelt file: " + ds.qmelt_file;
    } else{
        msg = "PROBLEM loading of Qmelt file: " + ds.qmelt_file;       
    } 
    std::cout << msg  << std::endl;
    logFLUXOSfile << msg + "\n";
    
    float tim = tmelts;
    return tim;
}