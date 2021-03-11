
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

#ifndef GLOBVARH_INCLUDED
#define GLOBVARH_INCLUDED

#include "jnlohmann/json.h"
using json = nlohmann::json;

class GlobVar
{
public:
  GlobVar() {

  }
  GlobVar(size_t MROWS, size_t MCOLS) {
    this->MROWS = MROWS;
    this->MCOLS = MCOLS;

    z= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    zb= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    h= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    ux= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    uy= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    qx= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    qy= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    qxf= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    qyf= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    us= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    dh= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    dqx = std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    dqy= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    //sbMROWS= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    //sbMCOLS= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    //cfri= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    ks= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    
//  f_1= mass flux per unit width
//  f_2= momentum flux per unit width in x-direction
//  f_3= momentum flux per unit width in y-direction
    fe_1= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    fe_2= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    fe_3= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    fn_1= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    fn_2= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    fn_3= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    twetimetracer= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS)); // connectivity-hours
    ldry= std::unique_ptr<arma::Mat<float>>( new  arma::fmat(MROWS,MCOLS));
    innerNeumannBCWeir= std::unique_ptr<arma::Mat<float>>( new  arma::fmat(MROWS,MCOLS));
        
    conc_SW= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    soil_mass= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    h0= std::unique_ptr<arma::Mat<double>>( new  arma::mat(MROWS,MCOLS));
    ldry_prev= std::unique_ptr<arma::Mat<float>>( new  arma::fmat(MROWS,MCOLS));
    
    meteo = std::unique_ptr<arma::Mat<float>>( new  arma::fmat(2000,2));
    inflow = std::unique_ptr<arma::Mat<float>>( new  arma::fmat(2000,2));
  }
    size_t NROWS,NCOLS;
    size_t MROWS,MCOLS,dxy,arbase,
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
        //sbMROWS,sbMCOLS,                                  // for calc of weight of water (bed slope term) (solver_wet)
        ks, //cfri                                  // Friction (Chezy model is not being used for now)
        fe_1,fe_2,fe_3,fn_1,fn_2,fn_3,twetimetracer,
        conc_SW,h0,soil_mass,basin_dem; 
    std::unique_ptr<arma::Mat<float>> ldry,innerNeumannBCWeir,meteo,inflow,ldry_prev;   
    double hdry,                                    //minimum water depth
        dtfl,tim,                                   // timestep for flow computation
        D_coef,soil_release_rate,soil_conc_bckgrd,qmelvtotal, qmelv_inc, SWEmax, SWEstd;
    
    std::string dem_file,meteo_file,inflow_file,sim_purp;
    unsigned long inflow_ncol, inflow_nrow;

    double NODATA_VALUE;

    json master_MODSET;
    
};

#endif