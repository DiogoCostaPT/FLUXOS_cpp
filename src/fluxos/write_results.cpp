

// Copyright 1992: Cornel Beffa and Roland Faeh
// Copyright 2013: Kashif Shaad and Diogo Costa
// Copyright 2019: Diogo Costa

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
#include <armadillo>
#include <memory> 
#include <chrono>
#include <ctime> 

#include "GlobVar.h"
#include "write_results.h"


bool write_results(
    GlobVar& ds, 
    int print_tag, 
    std::chrono::duration<double> elapsed_seconds)
{

    unsigned int icol,irow;
    int a = 0;
    double ux;
    
    std::string tprint = ds.output_folder + "/" + std::to_string(print_tag); 
    std::string filext(".txt");
    tprint += filext;

    arma::mat filedataR(ds.NROWS*ds.NCOLS,16); 
    
    for(icol=1;icol<=ds.NCOLS;icol++)
    {
        for(irow=1;irow<=ds.NROWS;irow++)
        {
            if (((*ds.zb)(irow,icol)!=0.0f) && ((*ds.h).at(irow,icol)>ds.h_min_print))
            {
                filedataR(a,0) = irow;
                filedataR(a,1) = icol;

                filedataR(a,2) = icol * (ds.dxy) + (ds.XLLCORNER);
                filedataR(a,3) = irow * (ds.dxy) + (ds.YLLCORNER);
                
                filedataR(a,4) = (*ds.z).at(irow,icol);
                filedataR(a,5) = (*ds.z).at(irow,icol) - (*ds.zb).at(irow,icol);
                filedataR(a,6) = (*ds.ux).at(irow,icol);
                filedataR(a,7) = (*ds.uy).at(irow,icol);
                filedataR(a,8) = (*ds.qx).at(irow,icol)*ds.dxy;
                filedataR(a,9) = (*ds.qy).at(irow,icol)*ds.dxy;
                filedataR(a,10) = (*ds.us).at(irow,icol); 
                // Only prints conc_SW[0]
                // if using openwq, then all the other concentrations will be in openwq
                filedataR(a,11) = (*ds.conc_SW)[0].at(irow,icol); // adesolver
                filedataR(a,12) = (*ds.soil_mass).at(irow,icol); // adesolver
                filedataR(a,13) = (*ds.fe_1).at(irow,icol);
                filedataR(a,14) = (*ds.fn_1).at(irow,icol);
                filedataR(a,15) = (*ds.twetimetracer).at(irow,icol);

                a = a + 1;

            }
        }
    }
   
    arma::mat filedata(std::max(0,a-1),15); 
    filedata = filedataR(arma::span(0,std::max(0,a-1)),arma::span(0,15));
    
    arma::field<std::string> header(filedata.n_cols);
    header(0) = "irow [-]";
    header(1) = "icol [-]";
    header(2) = "Xcoord";
    header(3) = "Ycoord";
    header(4) = "z [m]";
    header(5) = "h [m]";
    header(6) = "ux [m/s]";
    header(7) = "uy [m/s]";
    header(8) = "qx * dxy [m3/sec]";
    header(9) = "qy * dxy [m3/sec]";
    header(10) = "us [m/s]";
    header(11) = "conc_SW [mg/l]";
    header(12) = "soil_mass [g]";
    header(13) = "fe_1 [m2/s]";
    header(14) = "fn_1 [m2/s]";
    header(15) = "twetimetracer [sec]";

    bool outwritestatus =  filedata.save(arma::csv_name(tprint,header));
    return outwritestatus;
}
    