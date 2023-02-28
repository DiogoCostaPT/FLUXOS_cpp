
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

#include "GlobVar.h"
#include "initiate.h"
#include "common.h"

unsigned int initiation(
    GlobVar& ds,
    std::ofstream& logFLUXOSfile) {
    
    std::unique_ptr<double[]> zbs1(new double[ds.MROWS]);   
    double zbsw,zbnw,zbse,zbne,zbsum;
    unsigned int a,icol,irow,irow1,icol1;
    unsigned int NROWS1,NCOLS1;
    unsigned int timstart;
    NROWS1=ds.NROWS+1;
    NCOLS1=ds.NCOLS+1;

    // Create output file idf not existing
    check_mkdir(ds.output_folder);  

     // INTERPOLATE ELEVATIONS OF THE BOUNDARIES
    for(irow=0;irow<=NROWS1;irow++){
      zbs1[irow]=(*ds.zb).at(irow,1);
      (*ds.zb).at(irow,0) = (*ds.zb).at(irow,1);
      (*ds.zb).at(irow,NCOLS1) = (*ds.zb).at(irow,ds.NCOLS);
    }
    for(icol=0;icol<=NCOLS1;icol++){
      (*ds.zb).at(0,icol) = (*ds.zb).at(1,icol);
      (*ds.zb).at(NROWS1,icol) = (*ds.zb).at(ds.NROWS,icol); 
    }

    (*ds.zb).at(0,NCOLS1) = (*ds.zb).at(1,ds.NCOLS);
    (*ds.zb).at(NROWS1,NCOLS1) = (*ds.zb).at(ds.NROWS,ds.NCOLS);
    
    // INTERPOLATE CELL VERTICES TO CELL CENTER
    for(icol=1;icol<=ds.NCOLS;icol++){
      zbsw = zbs1[1];
      icol1  = icol+1;
      zbnw = (*ds.zb).at(1,icol1);
      zbs1[1]=zbnw;
      for(irow=1;irow<=ds.NROWS;irow++){
            irow1=irow+1;
            zbse =zbs1[irow1];
            zbne = (*ds.zb).at(irow1,icol1);
            a = 0;
            zbsum = 0;
            if (std::fabs(zbsw) != 99999){
                zbsum=zbsum + zbsw;
                a = a + 1;
            }
            if (std::fabs(zbse) != 99999){
                zbsum=zbsum + zbse;
                a = a + 1;
            }
            if (std::fabs(zbnw) != 99999){
                zbsum=zbsum + zbnw;
                a = a + 1;
            }
            if (std::fabs(zbne) != 99999){
                zbsum=zbsum + zbne;
                a = a + 1;
            }
            (*ds.zb).at(irow,icol)=(zbsum)/a;
            zbs1[irow1]=zbne;
            zbsw= zbse;
            zbnw= zbne;
      }
    }
        
    // INITIAL CONDITIONS
    for(icol=1;icol<=ds.NCOLS;icol++)
    {
        for(irow=1;irow<=ds.NROWS;irow++)
        {
            if(std::fabs((*ds.zb)(irow,icol))!=99999)
            {
            (*ds.h).at(irow,icol)=
                std::fmax((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0f);
            (*ds.z).at(irow,icol)=(*ds.zb).at(irow,icol)+(*ds.h).at(irow,icol);
            (*ds.qx).at(irow,icol)=(*ds.ux).at(irow,icol)*(*ds.h).at(irow,icol);
            (*ds.qy).at(irow,icol)=(*ds.uy).at(irow,icol)*(*ds.h).at(irow,icol);
            (*ds.soil_mass).at(irow,icol)  = ds.soil_conc_bckgrd;
            }
        }
    }
    
    ds.output_folder = ds.output_folder + "/";
    timstart = findLastStep(ds.output_folder.c_str()); // list the results files to get the last time step
    
    arma::mat filedata; 
    std::string init_file, msg;
    
    init_file = ds.output_folder + "/" + std::to_string(timstart) + ".txt";
    
    bool flstatus = filedata.load(init_file,arma::csv_ascii);

    if(flstatus == true) 
    {
        for(a=1;a<filedata.col(1).n_elem;a++)
        {
            irow = (int) (filedata(a,0));  
            icol = (int) (filedata(a,1)); 

            (*ds.z).at(irow,icol) = filedata(a,4);
            (*ds.h).at(irow,icol) = filedata(a,5); 
            (*ds.ux).at(irow,icol) = filedata(a,6);
            (*ds.uy).at(irow,icol) = filedata(a,7);
            (*ds.qx).at(irow,icol) = filedata(a,8);
            (*ds.qy).at(irow,icol) = filedata(a,9);
            (*ds.us).at(irow,icol) = filedata(a,10);
            (*ds.conc_SW).at(irow,icol) = filedata(a,11);
            (*ds.soil_mass).at(irow,icol) = filedata(a,12);
            (*ds.twetimetracer).at(irow,icol) = filedata(a,15);
            (*ds.ldry).at(irow,icol) = 0.0f;
        }
        msg = "Successful loading of initial conditions file: " + init_file;
    } else
    {
        msg = "NO INITIAL CONDITIONS FOUND: All variables set to zero";  

         for(icol=1;icol<=ds.NCOLS;icol++)
        {
            for(irow=1;irow<=ds.NROWS;irow++)
            {
                (*ds.h).at(irow,icol) = 0.0f;
                (*ds.z).at(irow,icol) = (*ds.zb).at(irow,icol);
                (*ds.ux).at(irow,icol) = 0.0f;
                (*ds.uy).at(irow,icol) = 0.0f;
                (*ds.qx).at(irow,icol) = 0.0f;
                (*ds.qy).at(irow,icol) = 0.0f;
                (*ds.us).at(irow,icol) = 0.0f;
                (*ds.ldry).at(irow,icol) = 1.0f;
            }
        }
    }
    std::cout << msg << std::endl;
    logFLUXOSfile << msg + "\n";
    
    // BOUNDARY VALUES (initial set up)
        for(icol=0;icol<=NCOLS1;icol++)
        {
          (*ds.zb).at(0,icol)=
            1.5*(*ds.zb).at(1,icol)-.5*(*ds.zb).at(2,icol);
          (*ds.zb).at(NROWS1,icol)=
            1.5*(*ds.zb).at(ds.NROWS,icol)-.5*(*ds.zb).at(ds.NROWS-1,icol);
          (*ds.z).at(0,icol)=
            1.5*(*ds.z).at(1,icol)-.5*(*ds.z).at(2,icol);
          (*ds.z).at(NROWS1,icol)=
            1.5*(*ds.z).at(ds.NROWS,icol)-.5*(*ds.z).at(ds.NROWS-1,icol);
          (*ds.h).at(0,icol)=
            std::fmax(0.0f,(*ds.z).at(0,icol)-(*ds.zb).at(0,icol));
          (*ds.h).at(NROWS1,icol)=
            std::fmax(0.0f,(*ds.z).at(NROWS1,icol)-(*ds.zb).at(NROWS1,icol));
          (*ds.qx).at(0,icol)=0.0f;
          (*ds.qy).at(NROWS1,icol)=0.0f;
          (*ds.qxf).at(0,icol)=0.0f;
        }
        for(irow=0;irow<=NROWS1;irow++)
        { 
          (*ds.zb).at(irow,0)=
            1.5*(*ds.zb).at(irow,1)-.5*(*ds.zb).at(irow,2);
          (*ds.zb).at(irow,NCOLS1)=
            1.5*(*ds.zb).at(irow,ds.NCOLS)-.5*(*ds.zb).at(irow,ds.NCOLS-1);
          (*ds.z).at(irow,0)=
            1.5*(*ds.z).at(irow,1)-.5*(*ds.z).at(irow,2);
          (*ds.z).at(irow,NCOLS1)=
            1.5*(*ds.z).at(irow,ds.NCOLS)-.5*(*ds.z).at(irow,ds.NCOLS-1);
          (*ds.h).at(irow,0)=
            std::fmax(0.0f,(*ds.z).at(irow,0)-(*ds.zb).at(irow,0));
          (*ds.h).at(irow,NCOLS1)=
            std::fmax(0.0f,(*ds.z).at(irow,NCOLS1)-(*ds.zb).at(irow,NCOLS1));
          (*ds.qx).at(irow,0)=0.0f;
          (*ds.qy).at(irow,NCOLS1)=0.0f;
          (*ds.qyf).at(irow,0)=0.0f;
        }
    
    return timstart;
}