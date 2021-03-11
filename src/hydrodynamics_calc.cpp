
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

#include <iostream>
#include <armadillo>
#include <memory> 

#include "GlobVar.h"
#include "hydrodynamics_calc.h"
#include "solver_drydomain.h"
#include "solver_wetdomain.h"

void hydrodynamics_calc(
    GlobVar& ds)
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
    float cell_neumann;

    dtl = ds.dtfl;
    
    // GET hp AND CHECK IF DRY OR WET
    for(icol=1;icol<=ds.NCOLS;icol++)
    {
        for(irow=1;irow<=ds.NROWS;irow++)
        {  
            hp=std::fmax(0.0f,(*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol));
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
        for(icol=1;icol<=ds.NCOLS;icol++)
        {
            for(irow=1;irow<=ds.NROWS;irow++)
            {  
                cell_neumann = (*ds.innerNeumannBCWeir).at(irow,icol);
                if((*ds.ldry).at(irow,icol) == 1.0f && cell_neumann == 0.0f) 
                {
                    solver_dry(ds,irow,icol);
                } else if (cell_neumann == 0.0f)
                {
                    solver_wet(ds,irow,icol);
                }
            }
        }
    }
    
    // CALCULATE TOTAL MASS AND MOMENTUM DERIVATIVE
    for(icol=1;icol<=ds.NCOLS;icol++)
    {
        for(irow=1;irow<=ds.NROWS;irow++)
        {  
             cell_neumann = (*ds.innerNeumannBCWeir).at(irow,icol);
             if (cell_neumann == 0.0f){
                (*ds.dh).at(irow,icol) = 
                    (((*ds.fe_1).at(irow-1,icol)-(*ds.fe_1).at(irow,icol))/ds.dxy 
                    +((*ds.fn_1).at(irow,icol-1)-(*ds.fn_1).at(irow,icol))/ds.dxy)*dtl;
                (*ds.dqx).at(irow,icol) = 
                    (((*ds.fe_2).at(irow-1,icol)-(*ds.fe_2).at(irow,icol))/ds.dxy 
                    +((*ds.fn_2).at(irow,icol-1)-(*ds.fn_2).at(irow,icol))/ds.dxy)*dtl;
                (*ds.dqy).at(irow,icol) = 
                    (((*ds.fe_3).at(irow-1,icol)-(*ds.fe_3).at(irow,icol))/ds.dxy 
                    +((*ds.fn_3).at(irow,icol-1)-(*ds.fn_3).at(irow,icol))/ds.dxy)*dtl;
                (*ds.qxf).at(irow,icol) = (*ds.fe_1).at(irow,icol)*dtl;
                (*ds.qyf).at(irow,icol) = (*ds.fn_1).at(irow,icol)*dtl;
             }else
             {
                (*ds.dh).at(irow,icol) = 0.0f;
                (*ds.dqx).at(irow,icol) = 0.0f;
                (*ds.dqy).at(irow,icol) = 0.0f;
                (*ds.qxf).at(irow,icol) = 0.0f;
                (*ds.qyf).at(irow,icol) = 0.0f;
             }  
        }
    }

    // CAL NEW VALUES
    for(icol=1;icol<=ds.NCOLS;icol++)
    {
        for(irow=1;irow<=ds.NROWS;irow++)
        {  
            (*ds.z).at(irow,icol)=(*ds.z).at(irow,icol)+(*ds.dh).at(irow,icol);
            hp=std::fmax(0.0f,(*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol));
            (*ds.h).at(irow,icol)=hp;
            cell_neumann = (*ds.innerNeumannBCWeir).at(irow,icol);
            
            if(hp<ds.hdry || cell_neumann == 1.0f) 
            {
                (*ds.qx).at(irow,icol)= 0.0f;
                (*ds.qy).at(irow,icol)= 0.0f;
                (*ds.us).at(irow,icol)= 0.0f;
                (*ds.ldry).at(irow,icol) = 1.0f;
            } else 
            {
                (*ds.qx).at(irow,icol)=
                    (*ds.qx).at(irow,icol)+(*ds.dqx).at(irow,icol);  // numerical flux at cell center
                (*ds.qy).at(irow,icol)=
                    (*ds.qy).at(irow,icol)+(*ds.dqy).at(irow,icol);  // numerical flux at cell center
                (*ds.qx).at(irow,icol)=
                    .1*(*ds.qxf).at(irow-1,icol)+.8*(*ds.qx).at(irow,icol)
                    +.1*(*ds.qxf).at(irow,icol); 
                (*ds.qy).at(irow,icol)=
                    .1*(*ds.qyf).at(irow,icol-1)+.8*(*ds.qy).at(irow,icol)
                    +.1*(*ds.qyf).at(irow,icol);
                (*ds.ldry).at(irow,icol) = 0.0f;          
            }  
        }
    } 
}
