
#ifndef GLOBVARH_INCLUDED
#define GLOBVARH_INCLUDED

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
    
    qmelt = std::unique_ptr<arma::Mat<float>>( new  arma::fmat(2000,2));
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
        conc_SW,h0,soil_mass; 
    std::unique_ptr<arma::Mat<float>> ldry,innerNeumannBCWeir,qmelt,ldry_prev;   
    double hdry,                                    //minimum water depth
        dtfl,tim,                                   // timestep for flow computation
        D_coef,soil_release_rate,soil_conc_bckgrd,qmelvtotal, qmelv_inc, SWEmax, SWEstd;
    
    std::string dem_file, qmelt_file,sim_purp;
    
};

#endif