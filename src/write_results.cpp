

#include <iostream>
#include <armadillo>
#include <memory> 
#include <chrono>
#include <ctime> 

#include "GlobVar.h"
#include "write_results.h"


bool write_results(GlobVar& ds, int print_tag, unsigned int print_step, std::chrono::duration<double> elapsed_seconds)
{

    unsigned int icol,irow;
    int a = 0;
    double ux;
    
    std::string tprint = "Results/" + std::to_string(print_tag); 
    std::string filext(".txt");
    tprint += filext;

    arma::mat filedataR(ds.NROWS*ds.NCOLS,14); 
    
    for(icol=1;icol<=ds.NCOLS;icol++)
    {
        for(irow=1;irow<=ds.NROWS;irow++)
        {
            if ((*ds.h).at(irow,icol)>0.0f)
            {
                filedataR(a,0) = irow;  
                filedataR(a,1) = icol; 
                filedataR(a,2) = (*ds.z).at(irow,icol); 
                filedataR(a,3) = (*ds.z).at(irow,icol) - (*ds.zb).at(irow,icol);
                filedataR(a,4) = (*ds.ux).at(irow,icol); 
                filedataR(a,5) = (*ds.uy).at(irow,icol); 
                filedataR(a,6) = (*ds.qx).at(irow,icol)*ds.dxy;
                filedataR(a,7) = (*ds.qy).at(irow,icol)*ds.dxy;
                filedataR(a,8) = (*ds.us).at(irow,icol); 
                filedataR(a,9) = (*ds.conc_SW).at(irow,icol); // adesolver
                filedataR(a,10) = (*ds.soil_mass).at(irow,icol); // adesolver
                filedataR(a,11) = (*ds.fe_1).at(irow,icol);
                filedataR(a,12) = (*ds.fn_1).at(irow,icol);
                filedataR(a,13) = (*ds.twetimetracer).at(irow,icol);
                a = a + 1;
            }
        }
    }
   
    arma::mat filedata(std::max(0,a-1),13); 
    filedata = filedataR(arma::span(0,std::max(0,a-1)),arma::span(0,13));
    
    arma::field<std::string> header(filedata.n_cols);
    header(0) = "irow [-]";
    header(1) = "icol [-]";
    header(2) = "z [m]";
    header(3) = "h [m]";
    header(4) = "ux [m/s]";
    header(5) = "uy [m/s]";
    header(6) = "qx * dxy [m3/sec]";
    header(7) = "qy * dxy [m3/sec]";
    header(8) = "us [m/s]";
    header(9) = "conc_SW [mg/l]";
    header(10) = "soil_mass [g]";
    header(11) = "fe_1 [m2/s]";
    header(12) = "fn_1 [m2/s]";
    header(13) = "twetimetracer [sec]";

    bool outwritestatus =  filedata.save(arma::csv_name(tprint,header));
    return outwritestatus;
}
    