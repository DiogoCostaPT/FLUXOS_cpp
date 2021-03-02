

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
#include "solver_drydomain.h"

void solver_dry(GlobVar& ds, unsigned int irow, unsigned int icol) {
    
    unsigned int iw,ie,is,in, NROWSl, NCOLSl;
    double fe1,fe2,fe3,fn1,fn2,fn3,zp,ze,zn,
           he,hn,qe,qp,rp,rn;
    double dze,hme,qme,dzn;
    double hmn,rmn,volrat;
    double ume,vmn,fe2p,fn3p;
    double hp;
    double zbp,zbe,zbn,zbpe,zbpn;
    double dtl;
    float ldw,ldp,lde,lds,ldn;
    double gaccl;
    
    is=icol-1;
    in=icol+1;
    iw=irow-1;
    ie=irow+1;
    
    ldw = (*ds.ldry).at(iw,icol);
    ldp = (*ds.ldry).at(irow,icol);
    lde = (*ds.ldry).at(ie,icol);
    lds = (*ds.ldry).at(irow,is);
    ldn = (*ds.ldry).at(irow,in);
    
    gaccl = ds.gacc;
    NCOLSl = ds.NCOLS;
    NROWSl = ds.NROWS;

    // CHECK IF ALL NEIGHBOUR CELLS ARE DRY
    if(ldw==1&&ldp==1&&lde==1&&lds==1&&ldn==1)
    {
        fe1=0.0f;
        fe2=0.0f;
        fe3=0.0f;
        fn1=0.0f;
        fn2=0.0f;
        fn3=0.0f;
        (*ds.dh).at(irow,icol)=0.0f;
        (*ds.dqx).at(irow,icol)=0.0f;
        (*ds.dqy).at(irow,icol)=0.0f;
        (*ds.qxf).at(irow,icol)=0.0f;
        (*ds.qyf).at(irow,icol)=0.0f;
        (*ds.qx).at(irow,icol)=0.0f;
        (*ds.qy).at(irow,icol)=0.0f;   
        return;
    }
    
    // CELL CENTER VALUES
    zbp = (*ds.zb).at(irow,icol);
    zbe = (*ds.zb).at(ie,icol);
    zbn = (*ds.zb).at(irow,in);
    zp=(*ds.z).at(irow,icol);
    ze=(*ds.z).at(ie,icol);
    zn=(*ds.z).at(irow,in);
    hp  = std::fmax(0.0f,(*ds.z).at(irow,icol)-(*ds.z).at(irow,icol));
    he=std::fmax(0.0f,ze-zbe);
    hn=std::fmax(0.0f,zn-zbn);
    qp=(*ds.qx).at(irow,icol);
    qe=(*ds.qx).at(ie,icol);
    rp=(*ds.qy).at(irow,icol);
    rn=(*ds.qy).at(irow,in);
       
    // CELL FACE VALUES  
    zbpe=.5*(zbe+zbp);
    zbpn=.5*(zbn+zbp);
    hme=.5*(hp+he);        
    qme=.5*(qp+qe);
    hmn=.5*(hp+hn);  
    rmn=.5*(rp+rn);
    dze=ze-zp;
    dzn=zn-zp;
   
    // TIMESTEP
    dtl = ds.dtfl;
    
    // INITIATION
    fe1=0.0f;
    fe2=0.0f;
    fe3=0.0f;
    fn1=0.0f;
    fn2=0.0f;
    fn3=0.0f;
    
    // CELLS WITH SOME DRY NEIGHBOURS
    if (lde==0.0f)
    {
        hme=std::fmax(0.0f,ze-zbpe);
        
        if(hme>ds.hdry) 
        {
            if(ze<=zp) 
            {
                fe2p=.5*gaccl*hme*hme;
                fe1=0.0f;
                fe2=fe2p;
                fe3=0.0f;
            }else 
            {
                dze=std::fmin(std::fabs(dze),hme);
                qme=0.296*dze*sqrt(gaccl*dze);      // Ritter solution
                hme=0.444*dze;                      // at cell side           
                volrat=ds.dxy*he/dtl;               // available volume rate per m of cell E [m2/s]
                qme=-std::fmin(qme,volrat);
                ume=qme/hme;                        // from cell center
                fe1=qme;
                fe2=qme*ume + .5*gaccl*hme*hme;
                fe3=0.0f;
            }
        }else 
        {
            fe1=0.0f;
            fe2=0.0f;
            fe3=0.0f;
        }
    }
    if (ldn==0.0f)  
    {   
        hmn=std::fmax(0.0f,zn-zbpn);    
        if(hmn>ds.hdry) 
        {
            if(zn<=zp) 
            {
                fn3p=.5*gaccl*hmn*hmn;
                fn1=0.0f;
                fn2=0.0f;
                fn3=fn3p;
            }else 
            {
            dzn=std::fmin(std::fabs(dzn),hmn);
            rmn=0.296*dzn*sqrt(gaccl*dzn);        // Ritter solution
            hmn=0.444*dzn;                       // at cell side           
            volrat=ds.dxy*hn/dtl;   // available volume rate per m of cell N [m2/s]
            rmn=-std::fmin(rmn,volrat);
            vmn=rmn/hmn;
            fn1=rmn;
            fn2=0.0f;
            fn3=rmn*vmn +.5*gaccl*hmn*hmn;
            }
        }else 
        {
            fn1=0.0f;
            fn2=0.0f;
            fn3=0.0f;
        }
    }
    
    // BOUNDARY CONDITIONS (WEIR DISCHARGE RATE)
    if (icol==1 || icol==NCOLSl)
    {
        fn1=std::fmin(volrat,sqrt(gaccl)*pow(std::fmax(hp,0.0f),1.5));
    }
    if (irow==1 || irow==NROWSl)
    {
        fe1=std::fmin(volrat,sqrt(gaccl)*pow(std::fmax(hp,0.0f),1.5));
    }
    
    // SAVE MASS AND MOMENTUM FLUXES
    (*ds.fn_1).at(irow,icol)=fn1;
    (*ds.fn_2).at(irow,icol)=fn2;
    (*ds.fn_3).at(irow,icol)=fn3;
    (*ds.fe_1).at(irow,icol)=fe1;
    (*ds.fe_2).at(irow,icol)=fe2;
    (*ds.fe_3).at(irow,icol)=fe3;   
} 
