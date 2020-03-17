

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

    arma::mat filedataR(ds.n_row*ds.n_col,15); 
    
    for(icol=1;icol<=ds.n_col;icol++)
    {
        for(irow=1;irow<=ds.n_row;irow++)
        {
            if ((*ds.h).at(irow,icol)>0.0f)
            {
                ux=sqrt((*ds.ux).at(irow,icol) * (*ds.ux).at(irow,icol) +
                        (*ds.uy).at(irow,icol) * (*ds.uy).at(irow,icol));
                filedataR(a,0) = irow;  
                filedataR(a,1) = icol; 
                filedataR(a,2) = (*ds.z).at(irow,icol); 
                filedataR(a,3) = (*ds.z).at(irow,icol) - (*ds.zb).at(irow,icol);
                filedataR(a,4) = (*ds.ux).at(irow,icol); 
                filedataR(a,5) = (*ds.uy).at(irow,icol); 
                filedataR(a,6) = (*ds.qx).at(irow,icol)*ds.dxy;
                filedataR(a,7) = (*ds.qy).at(irow,icol)*ds.dxy;
                filedataR(a,8) = ux; 
                filedataR(a,9) = (*ds.us).at(irow,icol); 
                filedataR(a,10) = (*ds.conc_SW).at(irow,icol); // adesolver
                filedataR(a,11) = (*ds.soil_mass).at(irow,icol); // adesolver
                filedataR(a,12) = (*ds.fe_1).at(irow,icol);
                filedataR(a,13) = (*ds.fn_1).at(irow,icol);
                filedataR(a,14) = (*ds.twetimetracer).at(irow,icol);
                a = a + 1;
            }
        }
    }
   
    arma::mat filedata(std::max(0,a-1),14); 
    filedata = filedataR(arma::span(0,std::max(0,a-1)),arma::span(0,14));
    
    bool outwritestatus =  filedata.save(tprint,arma::csv_ascii);
    return outwritestatus;
}
    