#include <armadillo>
#include <memory> 
#include <iostream>
#include <string>

#include "read_functions.h"

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
        if(i==1){ds.sim_purp = str;};
        if(i==2){ds.dem_file = pathfile + "/" + str;};
        if(i==3){ds.basin_file = pathfile + "/" + str;};
        if(i==4){ds.qmelt_file = pathfile + "/" + str;};
        if(i==5){(*print_step) = std::stoi(str);};
        if(i==6){(*ks_input) = std::stof(str);};  
        if(i==7){(ds.dxy) = std::stoi(str);};
        if(i==8){(ds.soil_release_rate) = std::stof(str);}; 
        if(i==9){(ds.soil_conc_bckgrd) = std::stof(str);};  
        if(i==10){(ds.SWEstd) = std::stof(str);}; 
        if(i==11){(ds.SWEmax) = std::stof(str);}; 
        
    }
    file.close();

    (ds.results_dir) = pathfile + "/Results";
    
    if(i==11){
        msg = "Successful loading of MODSET file: " + filename;
    } else{
        msg = "PROBLEM loading of MODSET file: " + filename;
    } 
     std::cout << msg  << std::endl;
     logFLUXOSfile << msg + "\n";
    
}

void read_geo(GlobVar& ds,double ks_input,std::ofstream& logFLUXOSfile)
{
    unsigned int icol,irow;  
    arma::mat filedata; 
    std::string msg;
    
    bool flstatus =  filedata.load(ds.dem_file,arma::raw_ascii);
 
   
    if(flstatus == true) {
        //for(a=0;a<filedata.col(1).n_elem;a++){
        for(icol=1;icol<=ds.n_col;icol++)
        {
            for(irow=1;irow<=ds.n_row;irow++)
            {   
                (*ds.zb).at(irow,icol) = std::abs(filedata(ds.n_row - irow,icol-1));  // For now -99999 is set as 99999 to act like a wall using abs
                (*ds.ks).at(irow,icol) = ks_input; 
            }
        }
        //}
        msg = "Successful loading of DEM file: " + ds.dem_file;
    } else{
        msg = "PROBLEM loading of DEM file: " + ds.dem_file;       
    } 
    std::cout << msg << std::endl;
    logFLUXOSfile << msg + "\n" ;
}

float read_load(GlobVar& ds,std::ofstream& logFLUXOSfile)
{
    unsigned int a; 
    unsigned int icol,irow;
    double tmelts,vmelt, tmelts_bef = 0.0f;
    std::string msg;
    
    arma::mat filedataB; 
    bool flstatusB =  filedataB.load(ds.basin_file,arma::raw_ascii);
    if(flstatusB == true) {
        for(icol=1;icol<=ds.n_col;icol++)
        {
            for(irow=1;irow<=ds.n_row;irow++)
            {  
                (*ds.basin_dem).at(irow,icol) = filedataB(ds.n_row - irow,icol-1);
            }  
        }
        msg = "Successful loading of Basin-DEM file: " + ds.basin_file;
    } else{
        msg = "PROBLEM loading of Basin-DEM file: " + ds.basin_file;       
    } 
    std::cout << msg  << std::endl;
    logFLUXOSfile << msg + "\n";
    
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