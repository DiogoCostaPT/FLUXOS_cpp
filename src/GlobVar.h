
#ifndef GLOBVARH_INCLUDED
#define GLOBVARH_INCLUDED

class GlobVar
{
public:
  GlobVar() {

  }
  GlobVar(size_t m_row, size_t m_col) {
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
    twetimetracer= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col)); // connectivity-hours
    ldry= std::unique_ptr<arma::Mat<float>>( new  arma::fmat(m_row,m_col));
        
    conc_SW= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    soil_mass= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    h0= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    ldry_prev= std::unique_ptr<arma::Mat<float>>( new  arma::fmat(m_row,m_col));
    
    basin_dem= std::unique_ptr<arma::Mat<double>>( new  arma::mat(m_row,m_col));
    qmelt = std::unique_ptr<arma::Mat<float>>( new  arma::fmat(2000,2));
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
        fe_1,fe_2,fe_3,fn_1,fn_2,fn_3,twetimetracer,
        conc_SW,h0,soil_mass,basin_dem; 
    std::unique_ptr<arma::Mat<float>> ldry,qmelt,ldry_prev;   
    double hdry,                                    //minimum water depth
        dtfl,tim,                                   // timestep for flow computation
        D_coef,soil_release_rate,soil_conc_bckgrd,qmelvtotal, qmelv_inc, SWEmax, SWEstd;
    
    std::string dem_file, basin_file, qmelt_file,sim_purp ;
    
};

#endif