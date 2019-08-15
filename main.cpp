
// Copyright 1992: Cornel Beffa and Roland Faeh
// Copyright 2013: Kashif Shaad and Diogo Costa
// Copyright 2019, Diogo Costa

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) an_col later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include<iostream>
#include<fstream>
#include<math.h>
#include<armadillo>
#include<string>
#include<memory> 
#include <chrono>
#include <ctime>  

#include <vector> 
#include <dirent.h>
#include <sys/types.h>

// read file names in Results directory
int findLastStep(const char *path) {

   struct dirent *entry;
   int i, timestart, filenum = 0, simnum;
   std::vector<char*> filenames; //stringvec filenames, filename_i;
   const char *filename_i;
   char *simnum_str_i;
   DIR *dir = opendir(path);
   
   if (dir != NULL) {
        while ((entry = readdir(dir)) != NULL) {
        filenames.push_back(entry->d_name); // storing the file names
        filenum = filenum + 1;
        }
   }
   closedir(dir);
   
   timestart = 0;
   for(i=2;i<filenum;i++){
       filename_i = filenames[i]; //.assign(filenames[i]); //strcpy(filename_i,(char *)(&filenames[i]));
        simnum_str_i = (char*) malloc(sizeof(filename_i)-2);
        strncpy (simnum_str_i, filename_i, sizeof(filename_i)-2);
        simnum = atoi(simnum_str_i);
        timestart = std::max(timestart,simnum);
        free(simnum_str_i);
   }
   
   //free(filename_i);
   free(entry);
   //free(dir);
   std::cout << "Start time (s): " << timestart << " (initial conditions available)" << std::endl;
   return timestart;
}

// get size of the domain
void get_domain_size(unsigned int *rown, unsigned int *coln )
{
    unsigned int numele;  
    arma::mat filedata; 
    bool flstatus =  filedata.load("dem_ersi_grid",arma::raw_ascii);
   
    *rown = 0;
    *coln = 0;
    
    if(flstatus == true) {
        *rown = filedata.col(1).n_elem;
        *coln = filedata.row(1).n_elem;
    } else{
        std::cout << "problem with loading 'dem_ersi_grid'" << std::endl;
    } 
}


class declavar
{
public:
  declavar() {

  }
  declavar(size_t m_row, size_t m_col) {
    this->m_row = m_row;
    this->m_col = m_col;

    z= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    zb= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    h= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    ux= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    uy= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    qx= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    qy= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    qxf= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    qyf= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    us= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    dh= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    dqx = std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    dqy= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    //sbm_row= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    //sbm_col= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    //cfri= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    ks= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    
//  f_1= mass flux per unit width
//  f_2= momentum flux per unit width in x-direction
//  f_3= momentum flux per unit width in y-direction
    fe_1= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    fe_2= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    fe_3= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    fn_1= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    fn_2= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    fn_3= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    ldry= std::unique_ptr<arma::Mat<float>>( new  arma::fmat(m_row,m_col));
    
    conc_SW= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    soil_mass= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    h0= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    ldry_prev= std::unique_ptr<arma::Mat<float>>( new  arma::fmat(m_row,m_col));
    
    basin_rowy= std::unique_ptr<arma::Mat<float>>( new  arma::fmat(210290,2));
    qmelt = std::unique_ptr<arma::Mat<float>>( new  arma::fmat(1633,2));
  }
    size_t n_row,n_col;
    size_t m_row,m_col,dxy,arbase, 
        ntim;                                       // maximum time step (seconds)
    float gacc = 9.80665,                           // gravitational acceleration 
        cfl,                                        // Courant condition
        cvdef, nuem;                                // for turbulence calc; nuem is molecular viscosity   
    std::unique_ptr<arma::Mat<double>> z,zb,h,
        ux,uy,                                        // velocities
        qx,qy,                                        // discharge at cell center: u*h [m3/s/m]
        qxf,qyf,                                      // discharges at cell face: u*h [m3/s/m]
        us,                                         // shear stress velocity 
        dh,dqx ,dqy,                                  // changes in h[irow][icol], p[irow][icol] and q[irow][icol]
        //sbm_row,sbm_col,                                  // for calc of weight of water (bed slope term) (solver_wet)
        ks, //cfri                                  // Friction (Chezy model is not being used for now)
        fe_1,fe_2,fe_3,fn_1,fn_2,fn_3,
        conc_SW,h0,soil_mass; 
    std::unique_ptr<arma::Mat<float>> ldry,basin_rowy,qmelt,ldry_prev;   
    double hdry,                                    //minimum water depth
        dtfl,tim,                                   // timestep for flow computation
        D_coef,soil_release_rate,soil_conc_bckgrd,qmelvtotal, qmelv_inc, SWEmax, SWEstd;
};


void read_geo(declavar& ds,double ks_input)
{
    unsigned int icol,irow;  
    arma::mat filedata; 
    bool flstatus =  filedata.load("dem_ersi_grid",arma::raw_ascii);
   
    if(flstatus == true) {
        //for(a=0;a<filedata.col(1).n_elem;a++){
        for(icol=1;icol<=ds.n_col;icol++)
        {
            for(irow=1;irow<=ds.n_row;irow++)
            {   
                (*ds.zb).at(irow,icol) = std::abs(filedata(ds.n_row - irow,icol-1));  
                (*ds.ks).at(irow,icol) = ks_input; 
            }
        }
        //}
    } else{
            std::cout << "problem with loading the DEM: file 'dem_ersi_grid'" << std::endl;
    } 
}

float read_load(declavar& ds)
{
    unsigned int a; 
    int icolb,irowb;
    double tmelts,vmelt, tmelts_bef = 0.0f;
    
    arma::mat filedataB; 
    bool flstatusB =  filedataB.load("Basin_Info.fluxos",arma::csv_ascii);
    if(flstatusB == true) {
        for(a=0;a<filedataB.col(0).n_elem;a++){
                       irowb = filedataB(a,0);  
            icolb = filedataB(a,1);  
            (*ds.basin_rowy).at(a,0) = irowb;  
            (*ds.basin_rowy).at(a,1) = icolb;  
            //printf("%f\n",(*ds.basin_rowy).at(a,0));
            //printf("%f\n",(*ds.basin_rowy).at(a,1));
        }
    } else{
            std::cout << "problem with loading 'Basin_Info.fluxos'" << std::endl;
    } 
    
    // reading qmelt 
    ds.qmelvtotal  = 0;
    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load("Qmelt_info.fluxos",arma::csv_ascii);
    if(flstatusQ == true) {
        for(a=0;a<filedataQ.col(1).n_elem;a++){
            tmelts = filedataQ(a,0);  // t melt seconds
            vmelt = filedataQ(a,1);  // value of melt
            (*ds.qmelt).at(a,0) = tmelts;  
            (*ds.qmelt).at(a,1) = vmelt;
            ds.qmelvtotal += vmelt /(1000.*3600.*24.) * (tmelts - tmelts_bef); 
            tmelts_bef = tmelts;
        }
    } else{
            std::cout << "problem with loading 'Qmelt_info.fluxos'" << std::endl;
    }
    float tim = tmelts;
    return tim;
}

unsigned int initiation(declavar& ds) {
    
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
            if (zbsw != 9999){
                zbsum=zbsum + zbsw;
                a = a + 1;
            }
            if (zbse != 9999){
                zbsum=zbsum + zbse;
                a = a + 1;
            }
            if (zbnw != 9999){
                zbsum=zbsum + zbnw;
                a = a + 1;
            }
            if (zbne != 9999){
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
            (*ds.h).at(irow,icol)=std::max((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0);
            (*ds.z).at(irow,icol)=(*ds.zb).at(irow,icol)+(*ds.h).at(irow,icol);
            (*ds.qx).at(irow,icol)=(*ds.ux).at(irow,icol)*(*ds.h).at(irow,icol);
            (*ds.qy).at(irow,icol)=(*ds.uy).at(irow,icol)*(*ds.h).at(irow,icol);
            (*ds.soil_mass).at(irow,icol)  = ds.soil_conc_bckgrd;
        }
    }
    
    timstart = findLastStep("Results/"); // list the results files to get the last time step
    
    arma::mat filedata; 
    bool flstatus = filedata.load("Results/" + std::to_string(timstart) + ".txt",arma::csv_ascii);

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
            (*ds.ldry).at(irow,icol) = 0.0f;
        }
    } else
    {
        std::cout << "No initial conditions (output files '*.txt' not found). All variables set to zero.'" << std::endl;
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

void solver_dry(declavar& ds, unsigned int irow, unsigned int icol) {
    
    unsigned int iw,ie,is,in, n_rowl, n_coll;
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
    n_coll = ds.n_col;
    n_rowl = ds.n_row;

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
    hp  = std::max(0.0,(*ds.z).at(irow,icol)-(*ds.z).at(irow,icol));
    he=std::max(0.0,ze-zbe);
    hn=std::max(0.0,zn-zbn);
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
        hme=std::max(0.0,ze-zbpe);
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
        hmn=std::max(0.0,zn-zbpn);    
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
    if (icol==1 || icol==n_coll)
    {
        fn1=std::min(volrat,sqrt(gaccl)*pow(std::fmax(hp,0.0f),1.5));
    }
    if (irow==1 || irow==n_rowl)
    {
        fe1=std::min(volrat,sqrt(gaccl)*pow(std::fmax(hp,0.0f),1.5));
    }
    
    // SAVE MASS AND MOMENTUM FLUXES
    (*ds.fn_1).at(irow,icol)=fn1;
    (*ds.fn_2).at(irow,icol)=fn2;
    (*ds.fn_3).at(irow,icol)=fn3;
    (*ds.fe_1).at(irow,icol)=fe1;
    (*ds.fe_2).at(irow,icol)=fe2;
    (*ds.fe_3).at(irow,icol)=fe3;   
} 

void solver_wet(declavar& ds, unsigned int irow, unsigned int icol){

    unsigned int iw,ie, is,in,inn, n_rowl, n_coll,dx,dy;
    double fe1,fe2,fe3,fn1,fn2,fn3,zw,zp,ze,zs,
        zn,znn,hw,he,hs,hn,hnn,fe2p,fn3p,qe,qp,qn,rw,rp,re,rn;
    double dze,dqe,dre,
        hme,qme,rme,dzn;
    double dqn,drn,hmn,hne,qmn,rmn,volrat,ume,vmn,volpot,cf;
    double cme,cmn,vme,umn,//dzbewr,dzbnsr,
        cnp,cne,cnn,up,un,us0,use0,une,vp,ve,vw,vwn,ven,txye,txyn,fe1c,fe2c,fe3c,fn1c,fn2c,fn3c,dzea,dhea;
    double cc,c1,c2,c3,c1a,c2a,c3a,a11,a12,a21,a22,a31,a32,a33,fe1r,fe2r,fe3r,dzna,dhna,a13,a23,fn1r,fn2r,fn3r;
    double hp ,hp0; 
    double zbw,zbp,zbe,zbs,zbn,zbnn,zbpe,zbpn; //zbpw,zbps
    double dtl, hdryl,gaccl,kspl, nueml, cvdefl; 
    float ldp,lde,ldn;
    bool lroe;

    n_rowl = ds.n_row;
    n_coll = ds.n_col;
    hdryl = ds.hdry;
    gaccl = ds.gacc;
    kspl = (*ds.ks).at(irow,icol);
    nueml = ds.nuem;
    cvdefl = ds.cvdef;
    
    is=icol-1;
    in=icol+1;
    inn=fmin(icol+2,n_coll+1);
    lroe = true;
    iw=irow-1;
    ie=irow+1;
   
    dx  = ds.dxy;
    dy  = ds.dxy;
    
    ldp = (*ds.ldry).at(irow,icol);
    lde = (*ds.ldry).at(ie,icol);
    ldn = (*ds.ldry).at(irow,in);

    // CELL CENTER VALUES
    zbw = (*ds.zb).at(iw,icol);
    zbp = (*ds.zb).at(irow,icol);
    zbe = (*ds.zb).at(ie,icol);
    zbs = (*ds.zb).at(irow,is);
    zbn = (*ds.zb).at(irow,in);
    zbnn= (*ds.zb).at(irow,inn);
    zw=(*ds.z).at(iw,icol);
    zp=(*ds.z).at(irow,icol);
    ze=(*ds.z).at(ie,icol);
    zs=(*ds.z).at(irow,is);
    zn=(*ds.z).at(irow,in);
    znn=(*ds.z).at(irow,inn);
    qp=(*ds.qx).at(irow,icol);
    qe=(*ds.qx).at(ie,icol);
    qn=(*ds.qx).at(irow,in);
    rw=(*ds.qy).at(iw,icol);
    rp=(*ds.qy).at(irow,icol);
    re=(*ds.qy).at(ie,icol);
    rn=(*ds.qy).at(irow,in);

    // zbpw=.5*(zbw+zbp);
    zbpe=.5*(zbe+zbp);
    // zbps=.5*(zbs+zbp);
    zbpn=.5*(zbn+zbp);
    hp  = std::max(0.0,(*ds.z).at(irow,icol)-zbp);
    hp0 = std::max(std::max(hdryl,hp),kspl);
    hw=std::max(0.0,zw-zbw);
    he=std::max(0.0,ze-zbe);
    hs=std::max(0.0,zs-zbs);
    hn=std::max(0.0,zn-zbn);
    hnn=std::max(0.0,znn-zbnn);

    // CELL FACE VALUES                    
    hme=.5*(hp+he);          
    qme=.5*(qp+qe);
    rme=.5*(rp+re);
    hmn=.5*(hp+hn);   
    qmn=.5*(qp+qn);
    rmn=.5*(rp+rn);

    dze=ze-zp;
    dqe=qe-qp;
    dre=re-rp;
    dzn=zn-zp;
    drn=rn-rp;
    dqn=qn-qp;
    
    // TIMESTEP
    dtl=ds.dtfl;

    // CELLS WITH SOME DRY NEIGHBOURS     
    if(lde==1.0f) 
    { 
        hme=std::max(0.0,zp-zbpe);
        he=0.0f;;
        qe=0.0f;
        re=0.0f;
        if(hme>hdryl) 
        {
             if(ze>=zp) {
                qme=0.0f;
                rme=0.0f;
            }else 
            {
                dze=fmin(fabs(dze),hme);
                qme=0.296*dze*sqrt(gaccl*dze);    // Ritter solution
                volrat=0.5*dx*hp/dtl;            // available volume rate per m of cell P [m2/s]
                qme=fmin(qme,volrat);
                rme=0.0f;
            }
        }else 
        {
            qme=0.0f;
            rme=0.0f;
            hme=0.0f;            
        }
        lroe= false;
    }
    if(ldn==1.0f) 
    {
        hmn=std::max(0.0,zp-zbpn);
        hn=0.0f;
        qn=0.0f;
        rn=0.0f;
        if(hmn>hdryl) 
        {
            if(zn>=zp) 
            {
                qmn=0.0f;
                rmn=0.0f;
            }else 
            {
                dzn=fmin(fabs(dzn),hmn);
                rmn=0.296*dzn*sqrt(gaccl*dzn);    // Ritter solution
                qmn=0.0f;
                volrat=0.5*dy*hp/dtl;            // available volume rate per m of cell P [m2/s]
                rmn=fmin(rmn,volrat);
            }
        }else 
        {
            qmn=0.0f;
            rmn=0.0f;
            hmn=0.0f;
        }
        lroe=false;
    }

    // CALC TURBULENT STRESS
    cme=sqrt(gaccl*hme);
    cmn=sqrt(gaccl*hmn);
    ume=qme/std::max(hme,hdryl);
    vme=rme/std::max(hme,hdryl);
    umn=qmn/std::max(hmn,hdryl);
    vmn=rmn/std::max(hmn,hdryl);
    
    cnp=cvdefl*(*ds.us).at(irow,icol)*hp+nueml;
    cne=cvdefl*(*ds.us).at(ie,icol)*he+nueml;
    cnn=cvdefl*(*ds.us).at(irow,in)*hn+nueml;
    hne=.5*(cnp+cne)*sqrt(hp*he);
    hnn=.5*(cnp+cnn)*sqrt(hp*hn);

    up=qp/hp0;
    un=qn/std::max(std::max(hn,hdryl),(*ds.ks).at(irow,in));
    us0=(*ds.qx).at(irow,is)/std::max(std::max(hs,hdryl),(*ds.ks).at(irow,is));
    use0=(*ds.qx).at(ie,is)/std::max(std::max((*ds.h).at(ie,is),hdryl),(*ds.ks).at(ie,is));
    une=(*ds.qx).at(ie,in)/std::max(std::max((*ds.h).at(ie,in),hdryl),(*ds.ks).at(ie,in));
    vp=rp/hp0;
    ve=re/std::max(std::max(he,hdryl),(*ds.ks).at(ie,icol));
    vw=rw/std::max(std::max(hw,hdryl),(*ds.ks).at(iw,icol));
    vwn=(*ds.qy).at(iw,in)/std::max(std::max((*ds.h).at(iw,in),hdryl),(*ds.ks).at(iw,in));
    ven=(*ds.qy).at(ie,in)/std::max(std::max((*ds.h).at(ie,in),hdryl),(*ds.ks).at(ie,in));

    txye=hne*((ve-vp)/fabs(dx)+.25*(un+une-us0-use0)/dy);
    txyn=hnn*((un-up)/fabs(dy)+.25*(ve+ven-vw-vwn)/dx);

    // CALC OF CONVECTION FLUXES
    fe1c=qme;
    fe2c=qme*ume;
    fe3c=qme*vme -txye;
    fn1c=rmn;
    fn2c=rmn*umn -txyn;
    fn3c=rmn*vmn;

    // ROE's DISSIPASSION
    if(lroe) 
    { 
        if(hme>ds.hdry) 
        {
            dzea=fabs(dze);
            if(dzea>0.5*hme) 
            {
                dhea=fabs(he-hp);
                dzea=fmin(dzea,dhea);
                dze=std::copysign(dzea,dze);
            }
            cc=.25/cme;
        }else 
        {
            cc=0.0f;
        }
        
        c1=ume;
        c2=ume+cme;
        c3=ume-cme;
        c1a=fabs(c1);
        c2a=fabs(c2);
        c3a=fabs(c3);
        a11=c2*c3a-c2a*c3;
        a12=c2a-c3a;
        a21=c2*c3*(c3a-c2a);
        a22=c2a*c2-c3a*c3;
        a31=vme*(c2*c3a-2.*cme*c1a-c2a*c3);
        a32=vme*(c2a-c3a);
        a33=2.*cme*c1a;

        fe1r=-(a11*dze+a12*dqe)*cc;
        fe2r=-(a21*dze+a22*dqe)*cc;
        fe3r=-(a31*dze+a32*dqe+a33*dre)*cc;

        if(ldp==0.0f&&hmn>ds.hdry) 
        {
                dzna=fabs(dzn);
                dhna=fabs(hn-hp);
                dzna=fmin(dzna,dhna);
                dzn=std::copysign(dzna,dzn);
                cc=.25/cmn;
        }else 
        {
            cc=0.0f;
        };
        c1=vmn;
        c2=vmn+cmn;
        c3=vmn-cmn;
        c1a=std::fabs(c1);
        c2a=std::fabs(c2);
        c3a=std::fabs(c3);
        a11=c2*c3a-c2a*c3;
        a13=c2a-c3a;
        a21=umn*(c2*c3a-2.*cmn*c1a-c2a*c3);
        a22=2.*cmn*c1a;
        a23=umn*(c2a-c3a);
        a31=c2*c3*(c3a-c2a);
        a33=c2a*c2-c3a*c3;

        fn1r=-(a11*dzn+a13*drn)*cc;
        fn2r=-(a21*dzn+a22*dqn+a23*drn)*cc;
        fn3r=-(a31*dzn+a33*drn)*cc;
    } else
    {
        fe1r=0.0f;
        fe2r=0.0f;
        fe3r=0.0f;
        fn1r=0.0f;
        fn2r=0.0f;
        fn3r=0.0f;
    }
    
    // PRESSURE AT CELL SIDE
    fe2p=.5*gaccl*(hme*hme);
    fn3p=.5*gaccl*(hmn*hmn); 

    // SUM OF ALL FLUXES
    fe1=fe1c+fe1r;
    fe2=fe2c+fe2r+fe2p;
    fe3=fe3c+fe3r;
    fn1=fn1c+fn1r;
    fn2=fn2c+fn2r;
    fn3=fn3c+fn3r+fn3p;
            
    // BOUNDARY CONDITIONS (WEIR DISCHARGE RATE) 
    if (icol==1 || icol==n_coll)
    {
        fn1=std::min(volrat,sqrt(gaccl)*pow(std::fmax(hp,0.0f),1.5));
    }
    if (irow==1 || irow==n_rowl)
    {
        fe1=std::min(volrat,sqrt(gaccl)*pow(std::fmax(hp,0.0f),1.5));
    }

    // CHECK MASS BALANCE (restrict outflow flux to available water)        
    volpot=ds.arbase*hp;    // volume in cell P [m3]
    volrat=volpot/dtl;      // max flux rate 
        
    if(volrat>0.0f) { // cell has water
        if(fe1>0.0f&&fn1>0.0f) 
        {
            if(fe1*dy+fn1*dx>volrat) {
                cf=fn1*dx/(fe1*dy+fn1*dx);
                fe1=(1.-cf)*volrat/dy;            // [m3/s]
                fn1=cf*volrat/dx;
            }
        }else if(fe1>0.0f) 
        {
            fe1=fmin(fe1*dy,volrat)/dy; 
        }else if(fn1>0.0f) 
        {
            fn1=fmin(fn1*dx,volrat)/dx;
        }
    }else // cell has no water 
    {
        fe1=0.0f;
        fn1=0.0f;
    }

    // SAVE MASS AND MOMENTUM FLUXES
    (*ds.fn_1).at(irow,icol)=fn1;
    (*ds.fn_2).at(irow,icol)=fn2;
    (*ds.fn_3).at(irow,icol)=fn3;
    (*ds.fe_1).at(irow,icol)=fe1;
    (*ds.fe_2).at(irow,icol)=fe2;
    (*ds.fe_3).at(irow,icol)=fe3;
    
} 


void flow_solver(declavar& ds)
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

// ADE solver
void adesolver(declavar& ds, int it)
{

    arma::mat qfcds(ds.m_row*ds.m_col,1);  //double qfcds(0:mx);
    arma::mat con_step(ds.m_row,ds.m_col);  //double qfcds(0:mx);
    double pfw,pfe,qfs,qfn,ntp, pfce, he,fp,fe, hne, pfde,area,areae,arean,hn,qxl,qyl,fw,
       fee,fs,fn,fnn,hnue,fem,hnn,qfcn,qfdn,fnm, cvolrate,cf,cbilan,dc,cvolpot,cvolrat,con, hnew;
    long unsigned int ix,iy,a;//!, printlim
    arma::mat cmaxr(ds.m_row,ds.m_col); //double  cmaxr(0:mx,0:my)
    arma::mat cminr(ds.m_row,ds.m_col); //cminr(0:mx,0:my);
    double dx,dy,dyn,hp,ie,iee,in, inn, is,iw;
    double nt =1 ; // eddy viscosity (m2/s) = 1,
    double sigc = 0.5;


    if(it>1) {
    // ADJUST CONCENTRATION TO NEW DEPTH
        for (a=1;a<=ds.n_col*ds.n_row;a++) {
            iy= ((a-1)/ds.n_row)+1;
            ix=a-ds.n_row*(iy-1);

            cmaxr(ix,iy)=std::max((*ds.conc_SW)(ix-1,iy),std::max((*ds.conc_SW)(ix+1,iy),std::max((*ds.conc_SW)(ix,iy-1),(*ds.conc_SW)(ix,iy+1))));
            cminr(ix,iy)=std::min((*ds.conc_SW)(ix-1,iy),std::min((*ds.conc_SW)(ix+1,iy),std::min((*ds.conc_SW)(ix,iy-1),(*ds.conc_SW)(ix,iy+1))));
            hnew=(*ds.h)(ix,iy);
            
            if((*ds.ldry)(ix,iy)==0 && (*ds.ldry_prev)(ix,iy)==0) 
            {
                (*ds.conc_SW)(ix,iy)=(*ds.conc_SW)(ix,iy)*(*ds.h0)(ix,iy)/hnew;
            } else if ((*ds.ldry)(ix,iy)==1) 
            {
                (*ds.conc_SW)(ix,iy) = 0.0f;
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
    dyn=ds.dxy; 

    // SPACE LOOP
    for (a=1;a<ds.n_col*ds.n_row;a++) {

        iy= ((a-1)/ds.n_row)+1;
        ix=a-ds.n_row*(iy-1);

        is=iy-1; 
        in=iy+1; 
        inn=std::min(iy+2,ds.n_col+1);
        iw=ix-1;
        ie=ix+1;
        iee=std::min(ix+2,ds.n_row+1);
                     
        //  BC 
        if (ix==1) {
            pfce=(*ds.conc_SW)(0,iy)*(*ds.fe_1)(0,iy)*dy;     // convective flux
            hp=std::max((*ds.h)(1,iy),ds.hdry);                  
            he=std::max((*ds.h)(2,iy),ds.hdry);
            fp=(*ds.conc_SW)(0,iy);
            fe=(*ds.conc_SW)(1,iy);
           
            hne=std::sqrt(hp*nt*he*nt)/sigc/std::abs(dx)*dy*ds.D_coef;
            pfde=0.;            // no diffusive flux over boundary
            pfe=pfce;  
        }  

        // CHECK IF THE DOMAIN IS DRY
        if((*ds.ldry)(ix,iy)==1){
            pfe=0.;
            qfcds(ix)=0.;
            (*ds.conc_SW)(ix,iy)=0.;            
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
        fw=(*ds.conc_SW)(iw,iy);
        fp=(*ds.conc_SW)(ix,iy);
        fe=(*ds.conc_SW)(ie,iy);
        fee=(*ds.conc_SW)(iee,iy);
        fs=(*ds.conc_SW)(ix,is);
        fn=(*ds.conc_SW)(ix,in);
        fnn=(*ds.conc_SW)(ix,inn);
            
        // FLUXES OVER WEST AND SOUTH FACES (from previous interaction)
        pfw=pfe; 
        qfs=qfcds(ix);


        // X-DIRECTION
        //// diffusive flux and mean concentration at east face
        if((*ds.ldry)(ie,iy)==0) {
            hnue=std::max(hp*nt*he*nt,.0001); 
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

        fem=std::max(0.,fem);

        if(ix==ds.n_row){  // if Boundary (overwrite the BC)
            fem=(*ds.conc_SW)(ds.n_row+1,iy);
        }

        //// advective flux - X-direction  - [m3/s]   
        pfce=qxl*fem*dy;  
        
        //// total flux = advective flux + diffusive
        pfe=pfce+pfde;      
        
        //// check available material if coming from the east cell
        if(pfe<0){ 
            if((*ds.ldry)(ie,iy)==0)    {
                cvolrate=-(fe*he)*areae/ds.dtfl; 
                pfe=std::max(pfe,cvolrate); //limit to available material
            }else {
                pfe=0.;
            }
        }             

        // Y-DIRECTION
        //// diffusion at the present time step Y-direction (pfde, where "d" refers to diffusion)
        if((*ds.ldry)(ix,in)==0)           {
            hnue=std::max(.0001,hp*ntp*hn*nt);
            hnn=std::sqrt(hnue)/sigc/dy*dx*ds.D_coef;              // [m3/s]
            qfdn=-hnn*(fn-fp);                    // diffusive flux
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

        fnm=std::max(0.,fnm);

        //// if Boundary (overwrite BC)
        if(iy==ds.n_col)    {
            fnm=(*ds.conc_SW)(ix,ds.n_col+1);
        }

        //// advective flux - X-direction  
        qfcn=qyl*fnm*dx; // [g/s]

        //// total flux
        qfn=qfcn+qfdn;
        
        //// check available material if coming from the north cell
        if(qfn<0)    {
            if((*ds.ldry)(ix,in)==0)    {     
                cvolrate=-(fn*hn)*arean/ds.dtfl; 
                qfn=std::max(qfn,cvolrate);   //limit to available material
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
                pfe=std::min(pfe,(cvolrat-qfn));
            } else if(qfn>0.){
                qfn=std::min(qfn,(cvolrat-pfe));         
            }
        }else { // bilance outflow with inflow
            if(pfe>=0.  &&  qfn<0.)  { //restrict pfe to bilan
                cbilan=cvolrat-qfn;                
                if(cbilan>0.) {
                    pfe=std::min(pfe,cbilan);
                }else{
                    pfe=0.;
                }
            } else if(pfe<0.  &&  qfn>=0.)  {
                cbilan=cvolrat-pfe;
                if(cbilan>0.) {
                    qfn=std::min(qfn,cbilan);
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
        con= (*ds.conc_SW)(ix,iy) +  dc/hp;

        con=std::fmin(cmaxr(ix,iy),con);
        con=std::fmax(cminr(ix,iy),con);
        (*ds.conc_SW)(ix,iy) = con;

        qfcds(ix)=qfn;  // convective+diffusive flux
  
    }
}

void wintra(declavar& ds)
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
                //deltam = (*ds.soil_mass).at(irow,icol) * (1-f) * ds.soil_release_rate/3600 * ds.dtfl; // mass release
                deltam = (*ds.soil_mass).at(irow,icol) * ds.soil_release_rate/3600 * ds.dtfl; // mass release
                (*ds.soil_mass).at(irow,icol) = (*ds.soil_mass).at(irow,icol) - deltam;
                (*ds.conc_SW).at(irow,icol) = (*ds.conc_SW).at(irow,icol) + deltam/(hp*ds.arbase);
            }
        }
    }
}

bool write_results(declavar& ds, int print_tag, unsigned int print_step, std::chrono::duration<double> elapsed_seconds)
{

    unsigned int icol,irow;
    int a = 0;
    double ux;
    
    std::string tprint = "Results/" + std::to_string(print_tag); 
    std::string filext(".txt");
    tprint += filext;

    arma::mat filedataR(ds.n_row*ds.n_col,14); 
    
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
                a = a + 1;
            }
        }
    }
   
    arma::mat filedata(std::max(0,a-1),14); 
    filedata = filedataR(arma::span(0,std::max(0,a-1)),arma::span(0,13));
    
    bool outwritestatus =  filedata.save(tprint,arma::csv_ascii);
    return outwritestatus;
}
    
int main(int argc, char** argv) 
{   
    unsigned int n_rowl, n_coll, it = 0;
    unsigned int a, irow, icol, print_step, print_next, qmelt_rowi, timstart;
    double c0,v0,u0,hp, hpall, qmelti,ntim_days,ks_input; 
    bool outwritestatus;
    std::chrono::duration<double> elapsed_seconds;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
        
    // Create/Open log file
    std::ofstream logFLUXOSfile ("fluxos_run.log");
   
     // Get the size of the domain (nrow and ncol)
    get_domain_size(&n_rowl, &n_coll);
     // Input the duration of the simulation
    
    std::cout << "FLUXOS"  << std::endl;
    logFLUXOSfile << "FLUXOS \n";
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::cout << "Simulation started... " << std::ctime(&start_time)  << std::endl;
    logFLUXOSfile << "Simulation started... " << std::ctime(&start_time);
    
    // Initiate variables on the heap
    declavar ds(n_rowl+2,n_coll+2); 
    
    // input/read data
    ds.cfl = 1; // Courant condition
    // ds.dxy = 3; // grid size (structure grid) - it will actually come from DEM
    ds.ntim = 0;// maximum time step (seconds)
    //kapa = -2.    // /  -2=1.Ord ; -1=2.Ord   // KOMISCH, DASS REAL/INTEGER ->schauen bei Rolands Dateien
    //betas = 2. // Chezy (parameter)
    //ksfirow = 0.2 // Chezy (rougness) -> NEEDs to be converted into a vector with data for all cells
    ds.cvdef = 0.07; // for turbulent stress calc
    ds.nuem = 1.2e-6; // molecular viscosity (for turbulent stress calc)
    //print_step = 3600; // in seconds
    // timstart = 558000; // start of the simulation
        
    
    // Request user input
    std::cout << "Print step (s) = ";
    std::cin >> print_step;
    logFLUXOSfile << "Print step (s) = " + std::to_string(print_step) + "\n";
  
    ds.n_row = ds.m_row - 2;
    ds.n_col = ds.m_col - 2;
    
    ds.D_coef = 0.01;
    
    std::cout << "Roughness height (m) = ";
    std::cin >> ks_input;
    logFLUXOSfile << "Roughness height (m) = " + std::to_string(ks_input) + "\n";
    
    std::cout << "Cell size (m) = ";
    std::cin >> ds.dxy;
    logFLUXOSfile << "Cell size (m) = " + std::to_string(ds.dxy) + "\n";
    
    ds.arbase = ds.dxy * ds.dxy;
    read_geo(ds,ks_input); // DEM
    ds.ntim = read_load(ds); // snowmelt load
    
    std::cout << "Simulation time (days) (Snowmelt input duration = " + std::to_string(ds.ntim/(3600*24)) + " days) = ";
    std::cin >> ntim_days;
    ds.ntim = ntim_days * 3600 * 24;
    logFLUXOSfile << "Simulation time (days) = " + std::to_string(ntim_days) + " (= " + std::to_string(ds.ntim) + " sec)";
   
    // Input the soil nutrient release rate
    std::cout << "Soil release rate (1/hour) = ";
    std::cin >> ds.soil_release_rate;
    logFLUXOSfile << "\nSoil release rate (1/hour) = " + std::to_string(ds.soil_release_rate);
    
    // Input the soil background concentration
    std::cout << "Soil initial background mass available for release to runoff (g) (0.txt points will be overwritten) = ";
    std::cin >> ds.soil_conc_bckgrd;
    logFLUXOSfile << "\nSoil initial background mass available for release to runoff (g) (0.txt points will be overwritten) = " + std::to_string(ds.soil_conc_bckgrd);
    
    std::cout << "SWE max (cm) = ";
    std::cin >> ds.SWEmax;
    logFLUXOSfile << "\nSWE max (cm) = " + std::to_string(ds.SWEmax);
    ds.SWEmax = ds.SWEmax/100;
    std::cout << "SWE std (cm) = ";
    std::cin >> ds.SWEstd;
    logFLUXOSfile << "\nSWE std (cm) = " + std::to_string(ds.SWEstd);
    ds.SWEstd = ds.SWEstd/100;
    
    timstart = initiation(ds);
    
    // INITIATION
    ds.hdry = (*ds.ks).at(1,1);  // temporary but basically saying that nothing will move until it reaches roughness height
        
    print_next = timstart;
    ds.tim = timstart;
    //write_results(ds,std::round(print_next));
    
    print_next = print_next + print_step;
        
    std::cout << "-----------------------------------------------\n" << std::endl;
    logFLUXOSfile << "\n-----------------------------------------------\n" << std::endl;
    
    // TIME LOOP
    while(ds.tim <= ds.ntim) 
    {              
        ds.dtfl=9.e10;
        hpall = 0.0f;
        
        // SPACE LOOP
        for(icol=1;icol<=ds.n_col;icol++)
        {
            for(irow=1;irow<=ds.n_row;irow++)
            {
                hp = (*ds.h).at(irow,icol);
                (*ds.h0)(irow,icol) = hp; // adesolver
                (*ds.ldry_prev).at(irow,icol) = (*ds.ldry).at(irow,icol); // adesolver
                        
                if(hp>ds.hdry)
                {
                    (*ds.ldry).at(irow,icol)=0.0f;
                    hp=std::fmax((*ds.h).at(irow,icol),ds.hdry);
                    hpall = std::fmax(hpall,(*ds.h).at(irow,icol));
                    c0=sqrt(ds.gacc*(*ds.h).at(irow,icol));
                    u0=std::fmax(.000001,fabs((*ds.qx).at(irow,icol)/hp));
                    v0=std::fmax(.000001,fabs((*ds.qy).at(irow,icol)/hp));
                    ds.dtfl=fmin(fmin(ds.cfl*ds.dxy/(u0+c0),ds.cfl*ds.dxy/(v0+c0)),ds.dtfl);
                }else 
                {
                    (*ds.ldry).at(irow,icol)=1.0f;
                }
                ds.dtfl=fmin(print_next - ds.tim, ds.dtfl);
            }
        }
                
        ds.tim = ds.tim + ds.dtfl;
          
       
        // Qmelt load
        for (a=0;a<=(*ds.qmelt).col(0).n_elem;a++){
            qmelt_rowi = a;
            if ((*ds.qmelt).at(a,0) > ds.tim){       
                break;
            }
        }
        
        qmelti = (*ds.qmelt).at(qmelt_rowi,1)/(1000.*3600.*24.)*ds.dtfl;
        ds.qmelv_inc += qmelti;
        for (a=0;a<(*ds.basin_rowy).col(1).n_elem;a++){
            irow = (*ds.basin_rowy).at(a,0);
            icol = (*ds.basin_rowy).at(a,1);
            hp = std::max((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0); // adesolver hp before adding snowmelt  
            (*ds.z).at(irow,icol) = (*ds.z).at(irow,icol) + qmelti;   
            (*ds.h)(irow,icol)=std::max((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0);
            if ((*ds.h)(irow,icol) <= ds.hdry)
            {
                (*ds.ldry).at(irow,icol)=0.0f;
        }
            (*ds.h0)(irow,icol) = (*ds.h)(irow,icol);
            if (hp!=0.)
            {          
                (*ds.conc_SW)(irow,icol)=((*ds.conc_SW)(irow,icol)*hp+qmelti*0)/((*ds.h)(irow,icol)); //adesolver (adjustment for snowmelt)       
            }
          }
                
        // FLOW SOLVERS
        if (hpall!=0) 
        {
            it++;
            flow_solver(ds);
            adesolver(ds, it);
            wintra(ds);
        }
        
        // PRINT RESULTS
        if (ds.tim>=print_next)
        {
             end = std::chrono::system_clock::now();
             elapsed_seconds = end-start;
             
              outwritestatus = write_results(ds,std::round(print_next),print_step,elapsed_seconds);
             
               
            if(outwritestatus == true) 
            {
                std::cout << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(print_step/60) << " || Time elapsed (min): " << elapsed_seconds.count()/60 << std::endl;
                logFLUXOSfile << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(print_step/60) << " || Time elapsed (min): " << std::to_string(elapsed_seconds.count()/60) + "\n";
                print_next = print_next + print_step;
                start = std::chrono::system_clock::now();

            } else
            {
                std::cout << "Problem when saving the results:" + print_next << std::endl;
                logFLUXOSfile << "Problem when saving the results:" + print_next;
                return 0;
            }
             
         }
    }
    
    // Simulation complete
    std::cout << "-----------------------------------------------" << std::endl;
    logFLUXOSfile << "\n-----------------------------------------------" << std::endl;
    std::cout << "Simulation complete (" << std::chrono::system_clock::now << ")"  << std::endl;
    logFLUXOSfile << "Simulation complete (" << std::chrono::system_clock::now;
    logFLUXOSfile.close(); 
      
    return 0;
}

