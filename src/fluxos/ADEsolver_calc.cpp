
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
#include "ADEsolver_calc.h"

// ADE solver
void adesolver_calc(
    GlobVar& ds, 
    int it,
    int ichem)
{

    arma::mat qfcds(ds.MROWS*ds.MCOLS,1);  //double qfcds(0:mx);
    arma::mat con_step(ds.MROWS,ds.MCOLS);  //double qfcds(0:mx);
    double pfw,pfe,qfs,qfn,ntp, pfce, he,fp,fe, hne, pfde,area,areae,arean,hn,qxl,qyl,fw,
       fee,fs,fn,fnn,hnue,fem,hnn,qfcn,qfdn,fnm, cvolrate,cf,cbilan,dc,cvolpot,cvolrat,con, hnew;
    long unsigned int ix,iy,a;//!, printlim
    arma::mat cmaxr(ds.MROWS,ds.MCOLS); //double  cmaxr(0:mx,0:my)
    arma::mat cminr(ds.MROWS,ds.MCOLS); //cminr(0:mx,0:my);
    double dx,dy,hp,ie,iee,in, inn, is,iw;
    double nt =1 ; // eddy viscosity (m2/s) = 1,
    double sigc = 0.5;


    if(it>1) {
    // ADJUST CONCENTRATION TO NEW DEPTH
        for (a=1;a<=ds.NCOLS*ds.NROWS;a++) {
            iy= ((a-1)/ds.NROWS)+1;
            ix=a-ds.NROWS*(iy-1);

            cmaxr(ix,iy)=std::fmax((*ds.conc_SW)[ichem](ix-1,iy),std::fmax((*ds.conc_SW)[ichem](ix+1,iy),std::fmax((*ds.conc_SW)[ichem](ix,iy-1),(*ds.conc_SW)[ichem](ix,iy+1))));
            cminr(ix,iy)=std::fmin((*ds.conc_SW)[ichem](ix-1,iy),std::fmin((*ds.conc_SW)[ichem](ix+1,iy),std::fmin((*ds.conc_SW)[ichem](ix,iy-1),(*ds.conc_SW)[ichem](ix,iy+1))));
            hnew=(*ds.h)(ix,iy);
            
            if((*ds.ldry)(ix,iy)==0 && (*ds.ldry_prev)(ix,iy)==0) 
            {
                (*ds.conc_SW)[ichem](ix,iy)=(*ds.conc_SW)[ichem](ix,iy)*(*ds.h0)(ix,iy)/hnew;
            } else if ((*ds.ldry)(ix,iy)==1) 
            {
                (*ds.conc_SW)[ichem](ix,iy) = 0.0f;
            }
        }
    }else
    {
        return;
    }
            
    //...    POLLUTION SOURCES
    ////$OMP PARALLEL
    //if (isqes2_on==1) call isqes2            // instantenous
    //if (csqes2_on==1) call csqes2          // continuous
    //if (gusqes2_on==1) call gusqes2        // gullies
    //if (gusqes2_on==1) call gusqes2_buildup
    //if (load_u_src>0) call usqes2                      // uniform
    ////$OMP END PARALLEL

    dx=ds.dxy;
    dy = ds.dxy;
    //dyn=ds.dxy; 

    // SPACE LOOP
    for (a=1;a<ds.NCOLS*ds.NROWS;a++) {

        iy= ((a-1)/ds.NROWS)+1;
        ix=a-ds.NROWS*(iy-1);

        is=iy-1; 
        in=iy+1; 
        inn=std::fmin(iy+2,ds.NCOLS+1);
        iw=ix-1;
        ie=ix+1;
        iee=std::fmin(ix+2,ds.NROWS+1);
                     
        //  BC 
        if (ix==1) {
            pfce=(*ds.conc_SW)[ichem](0,iy)*(*ds.fe_1)(0,iy)*dy;     // convective flux
            hp=std::fmax((*ds.h)(1,iy),ds.hdry);                  
            he=std::fmax((*ds.h)(2,iy),ds.hdry);
            fp=(*ds.conc_SW)[ichem](0,iy);
            fe=(*ds.conc_SW)[ichem](1,iy);
           
            hne=std::sqrt(hp*nt*he*nt)/sigc/std::fabs(dx)*dy*ds.D_coef;
            pfde=0.;            // no diffusive flux over boundary
            pfe=pfce;  
        }  

        // CHECK IF THE DOMAIN IS DRY
        if((*ds.ldry)(ix,iy)==1){
            pfe=0.;
            qfcds(ix)=0.;
            (*ds.conc_SW)[ichem](ix,iy)=0.;            
            continue; 
        };

        // INITIALIZATION 
        area=ds.arbase;
        areae=ds.arbase;
        arean=ds.arbase;
        ntp = nt; 
        hp=(*ds.h)(ix,iy); 
        he=(*ds.h)(ie,iy);
        hn=(*ds.h)(ix,in);
        qxl=(*ds.fe_1)(ix,iy);
        qyl=(*ds.fn_1)(ix,iy);
        fw=(*ds.conc_SW)[ichem](iw,iy);
        fp=(*ds.conc_SW)[ichem](ix,iy);
        fe=(*ds.conc_SW)[ichem](ie,iy);
        fee=(*ds.conc_SW)[ichem](iee,iy);
        fs=(*ds.conc_SW)[ichem](ix,is);
        fn=(*ds.conc_SW)[ichem](ix,in);
        fnn=(*ds.conc_SW)[ichem](ix,inn);
            
        // FLUXES OVER WEST AND SOUTH FACES (from previous interaction)
        pfw=pfe; 
        qfs=qfcds(ix);


        // X-DIRECTION
        //// diffusive flux and mean concentration at east face
        if((*ds.ldry)(ie,iy)==0) {
            hnue=std::fmax(hp*nt*he*nt,.0001); 
            hne=std::sqrt(hnue)/sigc/dx*dy*ds.D_coef; 
            pfde=-hne*(fe-fp);                     // diffusive flux

            if(qxl>0.0f){
                if ((*ds.ldry)(iw,iy)==0) {
                   fem=-.125*fw+.75*fp+.375*fe;
                }else {
                   fem=0.5*fp+0.5*fe;
                }             
            } else{
                if ((*ds.ldry)(iee,iy)==0) {
                  fem=.375*fp+.75*fe-.125*fee;
                }else {
                  fem=0.5*fp+0.5*fe;   
                }
            }
        }else {
            fem=0.;
            pfde=0.;
        }

        fem=std::fmax(0.,fem);

        if(ix==ds.NROWS){  // if Boundary (overwrite the BC)
            fem=(*ds.conc_SW)[ichem](ds.NROWS+1,iy);
        }

        //// advective flux - X-direction  - [m3/s]   
        pfce=qxl*fem*dy;  
        
        //// total flux = advective flux + diffusive
        pfe=pfce+pfde;      
        
        //// check available material if coming from the east cell
        if(pfe<0){ 
            if((*ds.ldry)(ie,iy)==0)    {
                cvolrate=-(fe*he)*areae/ds.dtfl; 
                pfe=std::fmax(pfe,cvolrate); //limit to available material
            }else {
                pfe=0.;
            }
        }             

        // Y-DIRECTION
        //// diffusion at the present time step Y-direction (pfde, where "d" refers to diffusion)
        if((*ds.ldry)(ix,in)==0)           {
            hnue=std::fmax(.0001,hp*ntp*hn*nt);
            hnn=std::sqrt(hnue)/sigc/dy*dx*ds.D_coef;              // [m3/s]
            qfdn=-hnn*(fn-fp);                                  // diffusive flux
            if(qyl>0.0f)       {
                    if((*ds.ldry)(ix,is)==0) {
                         fnm=-.125*fs+.75*fp+.375*fn; 
                    }else {
                        fnm=0.5*fp+.5*fn;    
                   }
            }else{
                   if ((*ds.ldry)(ix,inn)==0) {
                        fnm=.375*fp+.75*fn-.125*fnn;
                   }else {
                       fnm=.5*fp+.5*fn;
                  }
            }
        } else {
            fnm=0.;
            qfdn=0.;
        }

        fnm=std::fmax(0.,fnm);

        //// if Boundary (overwrite BC)
        if(iy==ds.NCOLS)    {
            fnm=(*ds.conc_SW)[ichem](ix,ds.NCOLS+1);
        }

        //// advective flux - X-direction  
        qfcn=qyl*fnm*dx; // [g/s]

        //// total flux
        qfn=qfcn+qfdn;
        
        //// check available material if coming from the north cell
        if(qfn<0)    {
            if((*ds.ldry)(ix,in)==0)    {     
                cvolrate=-(fn*hn)*arean/ds.dtfl; 
                qfn=std::fmax(qfn,cvolrate);   //limit to available material
            }else {
                qfn=0.;
            }
        }

        //// CHECK AVAILABLE MATERIAL - VOLUME RATE [m3/s] in actual cell
        cvolpot=(fp*hp)*area; // [m3]
        cvolrat=cvolpot/ds.dtfl + pfw+qfs; // inflow during actual time-step
      
        if (cvolrat>0.0f)  {                    // outflow is possible
            if(pfe>0.  &&  qfn>0.) {            // both outflow
                if (pfe+qfn > cvolrat) {        // limit outflow to volrat
                    cf=qfn/(pfe+qfn);
                    pfe= (1.-cf)*cvolrat;       // [m3/s]
                    qfn=cf*cvolrat;
                }
            } else if(pfe>0.){
                pfe=std::fmin(pfe,(cvolrat-qfn));
            } else if(qfn>0.){
                qfn=std::fmin(qfn,(cvolrat-pfe));         
            }
        }else { // bilance outflow with inflow
            if(pfe>=0.  &&  qfn<0.)  { //restrict pfe to bilan
                cbilan=cvolrat-qfn;                
                if(cbilan>0.) {
                    pfe=std::fmin(pfe,cbilan);
                }else{
                    pfe=0.;
                }
            } else if(pfe<0.  &&  qfn>=0.)  {
                cbilan=cvolrat-pfe;
                if(cbilan>0.) {
                    qfn=std::fmin(qfn,cbilan);
                }else{
                    qfn=0.;
                }
            } else if(pfe>=0.  &&  qfn>=0.)  {
                pfe=0.0f;
                qfn=0.0f;
            }        
        }                             

        // CALCULATE NEW CONCENTRATION
        dc=(pfw-pfe + qfs-qfn)*ds.dtfl/area;  // [m]  
        con= (*ds.conc_SW)[ichem](ix,iy) +  dc/hp;

        con=std::fmin(cmaxr(ix,iy),con);
        con=std::fmax(cminr(ix,iy),con);
        (*ds.conc_SW)[ichem](ix,iy) = con;

        qfcds(ix)=qfn;  // convective+diffusive flux
  
    }
}
