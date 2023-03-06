
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

#include "GlobVar.h"
#include "WINTRAsolver_calc.h"

void wintrasolver_calc(
    GlobVar& ds,
    int ichem)
{
    unsigned int icol,irow;
    double deltam,hp,zbp, f,frac, QVolstd;
    
    frac = ds.SWEmax/ds.SWEstd;
    QVolstd = ds.qmelvtotal/frac;
    
    f=tanh(1.26*(ds.qmelvtotal - ds.qmelv_inc)/QVolstd);
    
    
    for(icol=1;icol<=ds.NCOLS;icol++)
    {
        for(irow=1;irow<=ds.NROWS;irow++)
        {
            hp = (*ds.h).at(irow,icol);
            zbp =(*ds.zb).at(irow,icol);
            if(hp>ds.hdry & zbp != ds.NODATA_VALUE) 
            {       
                //deltam = (*ds.soil_mass).at(irow,icol) * (1-f) * ds.soil_release_rate/3600 * ds.dtfl; // exponential mass release
                //deltam = ds.soil_release_rate/3600 * ds.dtfl; // linear mass release 
                deltam = 
                    (*ds.soil_mass).at(irow,icol) * ds.soil_release_rate/3600 * ds.dtfl; // mass release
                (*ds.soil_mass).at(irow,icol) = 
                    (*ds.soil_mass).at(irow,icol) - deltam;
                (*ds.conc_SW)[ichem].at(irow,icol) = 
                    (*ds.conc_SW)[ichem].at(irow,icol) + deltam/(hp*ds.arbase);
            }
        }
    }
}