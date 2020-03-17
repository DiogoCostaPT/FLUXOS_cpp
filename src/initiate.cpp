
#include <armadillo>
#include <memory> 
#include <iostream>

#include "GlobVar.h"
#include "initiate.h"
#include "common.h"

unsigned int initiation(GlobVar& ds,std::ofstream& logFLUXOSfile) {
    
    std::unique_ptr<double[]> zbs1(new double[ds.m_row]);   
    double zbsw,zbnw,zbse,zbne,zbsum;
    unsigned int a,icol,irow,irow1,icol1;
    unsigned int n_row1,n_col1;
    unsigned int timstart;
    n_row1=ds.n_row+1;
    n_col1=ds.n_col+1;

     // INTERPOLATE ELEVATIONS OF THE BOUNDARIES
    for(irow=0;irow<=n_row1;irow++){
      zbs1[irow]=(*ds.zb).at(irow,1);
      (*ds.zb).at(irow,0) = (*ds.zb).at(irow,1);
      (*ds.zb).at(irow,n_col1) = (*ds.zb).at(irow,ds.n_col);
    }
    for(icol=0;icol<=n_col1;icol++){
      (*ds.zb).at(0,icol) = (*ds.zb).at(1,icol);
      (*ds.zb).at(n_row1,icol) = (*ds.zb).at(ds.n_row,icol); 
    }

    (*ds.zb).at(0,n_col1) = (*ds.zb).at(1,ds.n_col);
    (*ds.zb).at(n_row1,n_col1) = (*ds.zb).at(ds.n_row,ds.n_col);
    
    // INTERPOLATE CELL VERTICES TO CELL CENTER
    for(icol=1;icol<=ds.n_col;icol++){
      zbsw = zbs1[1];
      icol1  = icol+1;
      zbnw = (*ds.zb).at(1,icol1);
      zbs1[1]=zbnw;
      for(irow=1;irow<=ds.n_row;irow++){
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
    for(icol=1;icol<=ds.n_col;icol++)
    {
        for(irow=1;irow<=ds.n_row;irow++)
        {
            if(std::abs((*ds.basin_dem)(irow,icol))!=99999)
            {
            (*ds.h).at(irow,icol)=std::max((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0);
            (*ds.z).at(irow,icol)=(*ds.zb).at(irow,icol)+(*ds.h).at(irow,icol);
            (*ds.qx).at(irow,icol)=(*ds.ux).at(irow,icol)*(*ds.h).at(irow,icol);
            (*ds.qy).at(irow,icol)=(*ds.uy).at(irow,icol)*(*ds.h).at(irow,icol);
            (*ds.soil_mass).at(irow,icol)  = ds.soil_conc_bckgrd;
            }
        }
    }
    
    timstart = findLastStep("Results/"); // list the results files to get the last time step
    
    arma::mat filedata; 
    std::string init_file, msg;
    
    init_file = "Results/" + std::to_string(timstart) + ".txt";
    
    bool flstatus = filedata.load(init_file,arma::csv_ascii);

    if(flstatus == true) 
    {
        for(a=0;a<filedata.col(1).n_elem;a++)
        {
            irow = filedata(a,0);  
            icol = filedata(a,1);  
            (*ds.h).at(irow,icol) = filedata(a,3);
            (*ds.z).at(irow,icol) = (*ds.zb).at(irow,icol) + filedata(a,3);
            (*ds.ux).at(irow,icol) = filedata(a,4);
            (*ds.uy).at(irow,icol) = filedata(a,5);
            (*ds.qx).at(irow,icol) = filedata(a,6);
            (*ds.qy).at(irow,icol) = filedata(a,7);
            (*ds.us).at(irow,icol) = filedata(a,9);
            (*ds.conc_SW).at(irow,icol) = filedata(a,10);
            (*ds.soil_mass).at(irow,icol) = filedata(a,11);
            (*ds.twetimetracer).at(irow,icol) = filedata(a,14);
            (*ds.ldry).at(irow,icol) = 0.0f;
        }
        msg = "Successful loading of initial conditions file: " + init_file;
    } else
    {
        msg = "NO INITIAL CONDITIONS FOUND: All variables set to zero";  

         for(icol=1;icol<=ds.n_col;icol++)
        {
            for(irow=1;irow<=ds.n_row;irow++)
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
        for(icol=0;icol<=n_col1;icol++)
        {
          (*ds.zb).at(0,icol)=1.5*(*ds.zb).at(1,icol)-.5*(*ds.zb).at(2,icol);
          (*ds.zb).at(n_row1,icol)=1.5*(*ds.zb).at(ds.n_row,icol)-.5*(*ds.zb).at(ds.n_row-1,icol);
          (*ds.z).at(0,icol)=1.5*(*ds.z).at(1,icol)-.5*(*ds.z).at(2,icol);
          (*ds.z).at(n_row1,icol)=1.5*(*ds.z).at(ds.n_row,icol)-.5*(*ds.z).at(ds.n_row-1,icol);
          (*ds.h).at(0,icol)=std::max(0.0,(*ds.z).at(0,icol)-(*ds.zb).at(0,icol));
          (*ds.h).at(n_row1,icol)=std::max(0.0,(*ds.z).at(n_row1,icol)-(*ds.zb).at(n_row1,icol));
          (*ds.qx).at(0,icol)=0.0f;
          (*ds.qy).at(n_row1,icol)=0.0f;
          (*ds.qxf).at(0,icol)=0.0f;
        }
        for(irow=0;irow<=n_row1;irow++)
        { 
          (*ds.zb).at(irow,0)=1.5*(*ds.zb).at(irow,1)-.5*(*ds.zb).at(irow,2);
          (*ds.zb).at(irow,n_col1)=1.5*(*ds.zb).at(irow,ds.n_col)-.5*(*ds.zb).at(irow,ds.n_col-1);
          (*ds.z).at(irow,0)=1.5*(*ds.z).at(irow,1)-.5*(*ds.z).at(irow,2);
          (*ds.z).at(irow,n_col1)=1.5*(*ds.z).at(irow,ds.n_col)-.5*(*ds.z).at(irow,ds.n_col-1);
          (*ds.h).at(irow,0)=std::max(0.0,(*ds.z).at(irow,0)-(*ds.zb).at(irow,0));
          (*ds.h).at(irow,n_col1)=std::max(0.0,(*ds.z).at(irow,n_col1)-(*ds.zb).at(irow,n_col1));
          (*ds.qx).at(irow,0)=0.0f;
          (*ds.qy).at(irow,n_col1)=0.0f;
          (*ds.qyf).at(irow,0)=0.0f;
        }
    
    return timstart;
}