
#include <iostream>
#include <armadillo>
#include <memory> 

#include "GlobVar.h"
#include "hydrodynamics_calc.h"
#include "solver_drydomain.h"
#include "solver_wetdomain.h"

void hydrodynamics_calc(GlobVar& ds)
{
           
//-----------------------------------------------------------------------
// Solves shallow water equation for one time-step - Ritter solver
// Pressure term excluded from numerical flux
// Discretized as central difference. 
// Adjustment of source term
// MUSCL-approach with limiter in roe dissipation
// fw1= mass flux per unit width
// fw2= momentum flux per unit width in x-direction
// fw3= momentum flux per unit width in y-direction
//-----------------------------------------------------------------------

    unsigned int irow, icol;
    double hp, dtl;

    dtl = ds.dtfl;
    
    // GET hp AND CHECK IF DRY OR WET
    for(icol=1;icol<=ds.n_col;icol++)
    {
        for(irow=1;irow<=ds.n_row;irow++)
        {  
            hp=std::max(0.0,(*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol));
            (*ds.h).at(irow,icol) = hp;
            
            if(hp<=ds.hdry)
            {
              (*ds.qx).at(irow,icol)=0.0f;
              (*ds.qy).at(irow,icol)=0.0f;
              (*ds.us).at(irow,icol)=0.0f;
              (*ds.ldry).at(irow,icol) = 1.0f;;          
            } else
            {
              (*ds.ldry).at(irow,icol) = 0.0f;  
              (*ds.twetimetracer).at(irow,icol) += dtl/3600; 
            }
        }
    }
    
    // CALL FLOW SOLVERS (compute mass and momentum fluxes)
     #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for(icol=1;icol<=ds.n_col;icol++)
        {
            for(irow=1;irow<=ds.n_row;irow++)
            {  
                if((*ds.ldry).at(irow,icol) == 1.0f)
                {
                    solver_dry(ds,irow,icol);
                } else
                {
                    solver_wet(ds,irow,icol);
                }
            }
        }
    }
    
    // CALCULATE TOTAL MASS AND MOMENTUM DERIVATIVE
    for(icol=1;icol<=ds.n_col;icol++)
    {
        for(irow=1;irow<=ds.n_row;irow++)
        {  
            (*ds.dh).at(irow,icol)=(((*ds.fe_1).at(irow-1,icol)-(*ds.fe_1).at(irow,icol))/ds.dxy +((*ds.fn_1).at(irow,icol-1)-(*ds.fn_1).at(irow,icol))/ds.dxy)*dtl;
            (*ds.dqx).at(irow,icol)=(((*ds.fe_2).at(irow-1,icol)-(*ds.fe_2).at(irow,icol))/ds.dxy +((*ds.fn_2).at(irow,icol-1)-(*ds.fn_2).at(irow,icol))/ds.dxy)*dtl;
            (*ds.dqy).at(irow,icol)=(((*ds.fe_3).at(irow-1,icol)-(*ds.fe_3).at(irow,icol))/ds.dxy +((*ds.fn_3).at(irow,icol-1)-(*ds.fn_3).at(irow,icol))/ds.dxy)*dtl;
            (*ds.qxf).at(irow,icol)=(*ds.fe_1).at(irow,icol)*dtl;
            (*ds.qyf).at(irow,icol)=(*ds.fn_1).at(irow,icol)*dtl;
        }
    }

    // CAL NEW VALUES
    for(icol=1;icol<=ds.n_col;icol++)
    {
        for(irow=1;irow<=ds.n_row;irow++)
        {  
            (*ds.z).at(irow,icol)=(*ds.z).at(irow,icol)+(*ds.dh).at(irow,icol);
            hp=std::fmax(0.0f,(*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol));
            (*ds.h).at(irow,icol)=hp;
            
            if(hp<ds.hdry) 
            {
                (*ds.qx).at(irow,icol)= 0.0f;
                (*ds.qy).at(irow,icol)= 0.0f;
                (*ds.us).at(irow,icol)= 0.0f;
                (*ds.ldry).at(irow,icol) = 1.0f;
            } else 
            {
                (*ds.qx).at(irow,icol)=(*ds.qx).at(irow,icol)+(*ds.dqx).at(irow,icol);  // numerical flux at cell center
                (*ds.qy).at(irow,icol)=(*ds.qy).at(irow,icol)+(*ds.dqy).at(irow,icol);  // numerical flux at cell center
                (*ds.qx).at(irow,icol)=.1*(*ds.qxf).at(irow-1,icol)+.8*(*ds.qx).at(irow,icol)+.1*(*ds.qxf).at(irow,icol); 
                (*ds.qy).at(irow,icol)=.1*(*ds.qyf).at(irow,icol-1)+.8*(*ds.qy).at(irow,icol)+.1*(*ds.qyf).at(irow,icol);
                (*ds.ldry).at(irow,icol) = 0.0f;          
            }  
        }
    } 
}
