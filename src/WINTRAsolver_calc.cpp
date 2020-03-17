#include <iostream>
#include <armadillo>
#include <memory> 

#include "GlobVar.h"
#include "WINTRAsolver_calc.h"

void wintrasolver_calc(GlobVar& ds)
{
    unsigned int icol,irow;
    double deltam,hp,zbp, f,frac, QVolstd;
    
    frac = ds.SWEmax/ds.SWEstd;
    QVolstd = ds.qmelvtotal/frac;
    
    f=tanh(1.26*(ds.qmelvtotal - ds.qmelv_inc)/QVolstd);
    
    
    for(icol=1;icol<=ds.n_col;icol++)
    {
        for(irow=1;irow<=ds.n_row;irow++)
        {
            hp = (*ds.h).at(irow,icol);
            zbp =(*ds.zb).at(irow,icol);
            if(hp>ds.hdry && zbp != 9999) 
            {       
                //deltam = (*ds.soil_mass).at(irow,icol) * (1-f) * ds.soil_release_rate/3600 * ds.dtfl; // exponential mass release
                //deltam = ds.soil_release_rate/3600 * ds.dtfl; // linear mass release 
                deltam = (*ds.soil_mass).at(irow,icol) * ds.soil_release_rate/3600 * ds.dtfl; // mass release
                (*ds.soil_mass).at(irow,icol) = (*ds.soil_mass).at(irow,icol) - deltam;
                (*ds.conc_SW).at(irow,icol) = (*ds.conc_SW).at(irow,icol) + deltam/(hp*ds.arbase);
            }
        }
    }
}