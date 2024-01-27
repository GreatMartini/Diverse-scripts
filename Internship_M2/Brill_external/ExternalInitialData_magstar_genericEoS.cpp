
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <sys/stat.h>
#include <cstdlib>
#include <cstring>
#include "hdf5.h"

#include "ExternalInitialData.h"

using namespace std;
using namespace Lorene;

#define CUto10KM 1.476887220238973e-01
#define NBFM3_2_RHO_CGS 1.660000000000000e+14
#define DENSITY_CGS_2_SIM 1.619679366684996e-18
#define DENSITY_LOR_2_SIM 2.688667748697094e-04
#define MAG_FIELD_LOR_2_SIM 1.197536860265830e-07
#define NTABLES 6

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define CD4(paru, ip2, ip1, im1, im2, dx) ((((-paru[ip2]) + 8.0 * paru[ip1]) + ((-8.0 * paru[im1]) + paru[im2])) / (12.0 * dx))

#define CD6(paru, ip3, ip2, ip1, im1, im2, im3, dx) ((((paru[ip3]) - 9 * paru[ip2] + 45.0 * paru[ip1]) + ((-45.0 * paru[im1]) + 9 * paru[im2] - paru[im3])) / (60.0 * dx))
#define lessEq(a,b) ((fabs((a) - (b))/1.0E-15 > 10 ? false: (floor(fabs((a) - (b))/1.0E-15) < 1)) || (a)<(b))



double SAMRAI_External_Data::temp0, SAMRAI_External_Data::temp1;
double SAMRAI_External_Data::energy_shift;


double SAMRAI_External_Data::eos_rhomax, SAMRAI_External_Data::eos_rhomin;
double SAMRAI_External_Data::eos_tempmin, SAMRAI_External_Data::eos_tempmax;
double SAMRAI_External_Data::eos_yemin, SAMRAI_External_Data::eos_yemax;
double SAMRAI_External_Data::eos_epsmin, SAMRAI_External_Data::eos_epsmax;
double SAMRAI_External_Data::eos_pressmin, SAMRAI_External_Data::eos_pressmax;

double SAMRAI_External_Data::c2p_tempmin;
double SAMRAI_External_Data::c2p_tempmax;

int SAMRAI_External_Data::nrho;
int SAMRAI_External_Data::ntemp;
int SAMRAI_External_Data::nye;

double* SAMRAI_External_Data::alltables;
double* SAMRAI_External_Data::epstable;
double* SAMRAI_External_Data::presstable;
double* SAMRAI_External_Data::logrho;
double* SAMRAI_External_Data::logtemp;
double SAMRAI_External_Data::dlintemp,SAMRAI_External_Data::dlintempi;
double SAMRAI_External_Data::drholintempi;
double SAMRAI_External_Data::dlintempyei;
double SAMRAI_External_Data::drholintempyei;
double* SAMRAI_External_Data::yes;
double SAMRAI_External_Data::dtemp;
double SAMRAI_External_Data::dtempi;
double SAMRAI_External_Data::drho;
double SAMRAI_External_Data::drhoi;
double SAMRAI_External_Data::dye;
double SAMRAI_External_Data::dyei;
double SAMRAI_External_Data::drhotempi;
double SAMRAI_External_Data::drhoyei;
double SAMRAI_External_Data::dtempyei;
double SAMRAI_External_Data::drhotempyei;

/*double *eps,    double *tau,    double *rho,    double *D,   
  double *vx,     double *vy,     double *vz,   
  double *Sx,     double *Sy,     double *Sz,
  double *P,      double *sqcs, 
  double *Bx,     double *By,     double *Bz, 
  double *gxx,    double *gxy,    double *gxz,    double *gyy,  double *gyz,  double *gzz,
  double *Kxx,    double *Kxy,    double *Kxz,    double *Kyy,  double *Kyz,  double *Kzz, 
  double *Alpha,  double *chi,    double *trK,    double *theta,
  double *Gamh_x, double *Gamh_y, double *Gamh_z, 
  double *Betaux, double *Betauy, double *Betauz, 
  double Gamma,   double vacuum, double tov_angle, double tov_Asize, double tov_facp, 
  const double* dx, 
  double init_domain_x, double init_domain_y, double init_domain_z,   
  int init_patch_x,     int init_patch_y,     int init_patch_z, 
  const int d_ghost_width,
  int nx,               int ny,               int nz*/
void SAMRAI_External_Data::external_loadData(int vars, ...)
{

  va_list args;
  va_start(args, vars);

  int* eosType = va_arg(args, int*);
  double* initialYe = va_arg(args, double*);
  double* initialTemperature = va_arg(args, double*);
  double* minTableTemperature = va_arg(args, double*);
  double* minTableEnergy = va_arg(args, double*);
  double* energyshift = va_arg(args, double*);
  double* vacuum_ye_beta = va_arg(args, double*);
  double* vacuum_P_reset = va_arg(args, double*);
  double* vacuum_ye_reset = va_arg(args, double*);
  double* vacuum_temp_reset = va_arg(args, double*);
  double* vacuum_tau_reset = va_arg(args, double*);
  double* vacuum_rho_reset = va_arg(args, double*);

  double* vacuum_ye = va_arg(args, double*);
  double* vacuum_D = va_arg(args, double*);
  double* vacuum_tau = va_arg(args, double*);
  double* ye_maximum = va_arg(args, double*);
  double* temperature = va_arg(args, double *);
  double* Ye = va_arg(args, double *);
  double* DYf = va_arg(args, double *);
  double* eps = va_arg(args, double *);
  double* tau = va_arg(args, double *);
  double* rho = va_arg(args, double *);
  double* D = va_arg(args, double *);
  double* vx = va_arg(args, double *);
  double* vy = va_arg(args, double *);
  double* vz = va_arg(args, double *);
  double* Sx = va_arg(args, double *);
  double* Sy = va_arg(args, double *);
  double* Sz = va_arg(args, double *);
  double* P = va_arg(args, double *);
  double* sqcs = va_arg(args, double *);
  double* Bx = va_arg(args, double *);
  double* By = va_arg(args, double *);
  double* Bz = va_arg(args, double *);
  double* gxx = va_arg(args, double *);
  double* gxy = va_arg(args, double *);
  double* gxz = va_arg(args, double *);
  double* gyy = va_arg(args, double *);
  double* gyz = va_arg(args, double *);
  double* gzz = va_arg(args, double *);
  double* Kxx = va_arg(args, double *);
  double* Kxy = va_arg(args, double *);
  double* Kxz = va_arg(args, double *);
  double* Kyy = va_arg(args, double *);
  double* Kyz = va_arg(args, double *);
  double* Kzz = va_arg(args, double *);
  double* Alpha = va_arg(args, double *);
  double* chi = va_arg(args, double *);
  double* trK = va_arg(args, double *);
  double* theta = va_arg(args, double *);
  double* Gamh_x = va_arg(args, double *);
  double* Gamh_y = va_arg(args, double *);
  double* Gamh_z = va_arg(args, double *);
  double* Betaux = va_arg(args, double *);
  double* Betauy = va_arg(args, double *);
  double* Betauz = va_arg(args, double *);
  double* Gamma = va_arg(args, double*);
  double* vacuum_rho = va_arg(args, double*);
  double* tov_angle = va_arg(args, double*);
  double* tov_Asize = va_arg(args, double*);
  double* tov_facp = va_arg(args, double*);
  double* rho_0 = va_arg(args, double*);  
  double* rho_1 = va_arg(args, double*);  
  double* rho_2 = va_arg(args, double*);  
  double* a_0 = va_arg(args, double*);  
  double* a_1 = va_arg(args, double*);  
  double* a_2 = va_arg(args, double*);  
  double* a_3 = va_arg(args, double*);  
  double* K_0 = va_arg(args, double*);  
  double* K_1 = va_arg(args, double*);  
  double* K_2 = va_arg(args, double*);  
  double* K_3 = va_arg(args, double*);  
  double* gamma_0 = va_arg(args, double*);  
  double* gamma_1 = va_arg(args, double*);  
  double* gamma_2 = va_arg(args, double*);  
  double* gamma_3 = va_arg(args, double*);
  double* dx = va_arg(args, double *);
  double init_domain_x = va_arg(args, double);
  double init_domain_y = va_arg(args, double);
  double init_domain_z = va_arg(args, double);
  int init_patch_x = va_arg(args, int);
  int init_patch_y = va_arg(args, int);
  int init_patch_z = va_arg(args, int);
  int d_ghost_width = va_arg(args, int);
  int nx = va_arg(args, int);
  int ny = va_arg(args, int);
  int nz = va_arg(args, int);
  va_end(args);

  printf("$$ lorene_magstar: grid size nx=%d, ny=%d, nz=%d\n",
                      nx,ny,nz);

  const char * mag_data_name = "resu.d";

  struct stat buffer;   
  if (stat ("resu.d", &buffer) != 0) {
      std::cout <<"Lorene exception: File resu.d does not exist."<< std::endl;
      std::exit(-1);
  }

  cerr << "opening resu file" << endl;
  FILE * mag_data_file = fopen(mag_data_name, "r");
  cerr << "done opening resu file" << endl;
  if (mag_data_file == NULL) printf("file handle problem.  The data file %s cannot be found.\n",mag_data_name);
  

  Mg3d spectral_grid(mag_data_file);
  Map_et mapping(spectral_grid, mag_data_file);

  cerr << "calling read eos" << endl;
  Eos * p_eos = Eos::eos_from_file(mag_data_file);
  cerr << "done calling eos" << endl;

  Et_rot_mag star(mapping, *p_eos, mag_data_file);
  star.equation_of_state();
  star.update_metric();
  star.extrinsic_curvature();
  star.hydro_euler();
  
  Scalar sp_density = star.get_nbar()() ;
  Scalar sp_pressure = star.get_press()() ;

  Cmp sp_v_x = star.get_u_euler()(0) ; //contravariant representation!!
  Cmp sp_v_y = star.get_u_euler()(1) ;
  Cmp sp_v_z = star.get_u_euler()(2) ;

  const Map& mp = star.get_mp();

  const Base_vect_spher& bspher = mp.get_bvect_spher() ;

  Sym_tensor gij(mp, COV, bspher) ;
  gij.set_etat_zero() ;
  gij.set(1,1) = star.get_a_car()() ;
  gij.set(2,2) = star.get_a_car()() ;
  gij.set(3,3) = star.get_b_car()() ;
  Metric gam(gij) ; //the 3-metric in quasi-isotropic coordinates

  Scalar fac = sqrt(star.get_a_car()()) ;//to transform B^(i) into B^i
  fac.std_spectral_base() ;
  Vector Bmag(mp, CON, bspher) ;
  Bmag.set(1) = Scalar(star.Magn()(0)) / fac ;
  Bmag.set(1).dec_dzpuis(2) ;
  Bmag.set(2) = Scalar(star.Magn()(1)) / fac ;
  Bmag.set(2).dec_dzpuis(2) ;
  Bmag.set(3) = 0 ;
  Tbl maxB = 0.5*(max(abs(Bmag(1))) + max(abs(Bmag(2)))) ;
#ifndef __clang__
  cout << maxB << endl ;
  cout << "div(B) / max(B) in each domain :  " <<
    max(abs(Bmag.divergence(gam))) / maxB << endl;
#endif
  //the divergence using the 3-metric associated covariant derivative

  Bmag.change_triad(mp.get_bvect_cart()) ;
  Scalar sp_mag_x = Bmag(1);//contravariant representation!!
  Scalar sp_mag_y = Bmag(2);
  Scalar sp_mag_z = Bmag(3);

  Cmp sp_lapse = star.get_nnn()() ;//the lapse
  Cmp sp_shift_x = star.get_shift()(0) ; //contravariant representation!!
  Cmp sp_shift_y = star.get_shift()(1) ;
  Cmp sp_shift_z = star.get_shift()(2) ;


  //the spatial metric g_{ij}
  Tenseur_sym sp_gamma(mapping, 2, COV, mapping.get_bvect_spher()) ;
  sp_gamma.set_etat_qcq() ;

  for (int i=0; i<3; i++) {
    for (int j=0; j<i; j++)
      sp_gamma.set(i,j) = 0 ;
    if (i != 2) sp_gamma.set(i,i) = star.get_a_car()() ;
    else sp_gamma.set(2,2) = star.get_bbb()()*star.get_bbb()() ;
  }

  sp_gamma.change_triad(mapping.get_bvect_cart()) ;
  Cmp sp_gam_xx = sp_gamma(0,0) ;
  Cmp sp_gam_xy = sp_gamma(0,1) ;
  Cmp sp_gam_xz = sp_gamma(0,2) ;
  Cmp sp_gam_yy = sp_gamma(1,1) ;
  Cmp sp_gam_yz = sp_gamma(1,2) ;
  Cmp sp_gam_zz = sp_gamma(2,2) ;

  //the extrinsic curvature
  Tenseur_sym sp_kij = star.get_bbb()*star.get_bbb()*star.get_tkij() ;
  for (int i=0; i<3; i++)
    for (int j=0; j<=i; j++) {
      sp_kij.set(i,j).dec2_dzpuis() ;
    }

  Cmp sp_kij_xx = sp_kij(0,0) ; //covariant form...
  Cmp sp_kij_xy = sp_kij(0,1) ;
  Cmp sp_kij_xz = sp_kij(0,2) ;
  Cmp sp_kij_yy = sp_kij(1,1) ;
  Cmp sp_kij_yz = sp_kij(1,2) ;
  Cmp sp_kij_zz = sp_kij(2,2) ;

  //assign values to cactus variables on a cartesian grid
  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        int i3D = i + nx*(j  + ny*k);

        // convert position from cactus units to Lorene units (10Km)
        double x_fd = (init_domain_x + ((i + init_patch_x) - d_ghost_width) * dx[0])*CUto10KM;
        double y_fd = (init_domain_y + ((j + init_patch_y) - d_ghost_width) * dx[1])*CUto10KM;
        double z_fd = (init_domain_z + ((k + init_patch_z) - d_ghost_width) * dx[2])*CUto10KM;

        double r_sp = sqrt(x_fd*x_fd + y_fd*y_fd + z_fd*z_fd) ;
        double theta_sp = acos(z_fd / r_sp) ;

        double phi_sp = atan(y_fd/x_fd); //this goes from -pi/2 to +pi/2

        if (x_fd<0.0) phi_sp = 2.0*acos(0.0)+phi_sp; //pi+phi

        if (r_sp==0) theta_sp=0;
        if ((x_fd*x_fd+y_fd*y_fd)==0) phi_sp=0;

        // the lapse
        Alpha[i3D]   = sp_lapse.val_point(r_sp, theta_sp, phi_sp);

        // and the shift beta^i
        Betaux[i3D] = -sp_shift_x.val_point(r_sp, theta_sp, phi_sp);
        Betauy[i3D] = -sp_shift_y.val_point(r_sp, theta_sp, phi_sp);
        Betauz[i3D] = -sp_shift_z.val_point(r_sp, theta_sp, phi_sp);
 
      /* En had no se hace, pero se pensaba que si  
        Betaux[i3D] = 0;
        Betauy[i3D] = 0;
        Betauz[i3D] = 0;*/
 
        double numdensity = sp_density.val_point(r_sp, theta_sp, phi_sp) ;
        rho[i3D]   = NBFM3_2_RHO_CGS * DENSITY_CGS_2_SIM * numdensity;
        P[i3D]     = sp_pressure.val_point(r_sp, theta_sp, phi_sp);
        // Covert the various quantities to HAD units
        P[i3D]    *= DENSITY_LOR_2_SIM;

        vx[i3D]    = sp_v_x.val_point(r_sp, theta_sp, phi_sp);
        vy[i3D]    = sp_v_y.val_point(r_sp, theta_sp, phi_sp);
        vz[i3D]    = sp_v_z.val_point(r_sp, theta_sp, phi_sp);

        if (*eosType == 0) {
          
          if (P[i3D]   < (*vacuum_rho) * (*vacuum_rho)) P[i3D] = (*vacuum_rho) * (*vacuum_rho);
          if (rho[i3D] < *vacuum_rho) {
            rho[i3D] = *vacuum_rho;
            vx[i3D] = 0.0;
            vy[i3D] = 0.0;
            vz[i3D] = 0.0;
          } 

          double p_cold, eps_cold, selected_gamma;
          if (rho[i3D] > (*rho_2)) {
            selected_gamma = (*gamma_3);
            p_cold = (*K_3) * (pow(rho[i3D], (*gamma_3)));
            eps_cold = (*a_3) + (*K_3) / ((*gamma_3) - 1.0) * (pow(rho[i3D], ((*gamma_3) - 1.0)));
          } else {
            if (!((*rho_1) > rho[i3D])) {
              selected_gamma = (*gamma_2);
              p_cold = (*K_2) * (pow(rho[i3D], (*gamma_2)));
              eps_cold = (*a_2) + (*K_2) / ((*gamma_2) - 1.0) * (pow(rho[i3D], ((*gamma_2) - 1.0)));
            } else {
              if (!((*rho_0) > rho[i3D])) {
                selected_gamma = (*gamma_1);
                p_cold = (*K_1) * (pow(rho[i3D], (*gamma_1)));
                eps_cold = (*a_1) + (*K_1) / ((*gamma_1) - 1.0) * (pow(rho[i3D], ((*gamma_1) - 1.0)));
              } else {
                selected_gamma = (*gamma_0);
                p_cold = (*K_0) * (pow(rho[i3D], (*gamma_0)));
                eps_cold = (*a_0) + (*K_0) / ((*gamma_0) - 1.0) * (pow(rho[i3D], ((*gamma_0) - 1.0)));
              }
            }
          }

          eps[i3D]   = eps_cold + (P[i3D] - p_cold) / (((*Gamma) - 1.0) * rho[i3D]);

          double h = rho[i3D] * (1.0 + eps[i3D]) + P[i3D];
          sqcs[i3D] = fabs(((*Gamma) * P[i3D] + (selected_gamma - (*Gamma)) * p_cold) / h);
        }
        // B^j
        Bx[i3D]    = sp_mag_x.val_point(r_sp, theta_sp, phi_sp);
        By[i3D]    = sp_mag_y.val_point(r_sp, theta_sp, phi_sp);
        Bz[i3D]    = sp_mag_z.val_point(r_sp, theta_sp, phi_sp);

        // g_ij
        gxx[i3D]   = sp_gam_xx.val_point(r_sp, theta_sp, phi_sp);
        gxy[i3D]   = sp_gam_xy.val_point(r_sp, theta_sp, phi_sp);
        gxz[i3D]   = sp_gam_xz.val_point(r_sp, theta_sp, phi_sp);
        gyy[i3D]   = sp_gam_yy.val_point(r_sp, theta_sp, phi_sp);
        gyz[i3D]   = sp_gam_yz.val_point(r_sp, theta_sp, phi_sp);
        gzz[i3D]   = sp_gam_zz.val_point(r_sp, theta_sp, phi_sp);

        // the extrinsic curvature K_ij
        Kxx[i3D]  = sp_kij_xx.val_point(r_sp, theta_sp, phi_sp);
        Kxy[i3D]  = sp_kij_xy.val_point(r_sp, theta_sp, phi_sp);
        Kxz[i3D]  = sp_kij_xz.val_point(r_sp, theta_sp, phi_sp);
        Kyy[i3D]  = sp_kij_yy.val_point(r_sp, theta_sp, phi_sp);
        Kyz[i3D]  = sp_kij_yz.val_point(r_sp, theta_sp, phi_sp);
        Kzz[i3D]  = sp_kij_zz.val_point(r_sp, theta_sp, phi_sp);

        // eps doesn't need to be converted, same units in Cactus and Lorene
        // [K_ij] = 1/[L]
        Kxx[i3D]  *= CUto10KM;
        Kxy[i3D]  *= CUto10KM;
        Kxz[i3D]  *= CUto10KM;
        Kyy[i3D]  *= CUto10KM;
        Kyz[i3D]  *= CUto10KM;
        Kzz[i3D]  *= CUto10KM;
        // v^j don't need to be converted, they are already in c unit
        // (c is the speed of light)

        Bx[i3D] *= MAG_FIELD_LOR_2_SIM;
        By[i3D] *= MAG_FIELD_LOR_2_SIM;
        Bz[i3D] *= MAG_FIELD_LOR_2_SIM;

        theta[i3D] = 0;
      }
    }
  }

  if (*eosType == 1) {
    double xye, xeps, xcsq, xp;
    std:string filename = "InitialData.h5";
    nuc_eos_C_ReadTable(const_cast<char*>(filename.c_str()));
    //Setting Vacuums
    double xtemp = eos_tempmin;
    xye   = *vacuum_ye_beta;
    double xrho  = *vacuum_rho;

    energy_shift = *energyshift;


    nuc_eos_lookup(xp, xeps, xtemp, xcsq, xrho, xye);
    *vacuum_D = xrho;
    *vacuum_tau = xrho*xeps;
    *vacuum_ye = eos_yemin;


    std::cout<<"****************ExternalInitialData********************"<<std::setprecision(15)<<std::endl;
    std::cout<<"These parameters have to be copied to the parameter file for restarting"<<std::endl;
    std::cout<<"vacuum_tau = "<<*vacuum_tau<<std::endl;
    std::cout<<"vacuum_ye = "<<*vacuum_ye<<std::endl;
    std::cout<<"vacuum_ye_beta = "<<*vacuum_ye_beta<<std::endl;
    std::cout<<"vacuum_D = "<<*vacuum_D<<std::endl;

    xrho = 1.2 * xrho;
    xtemp = 1.2 * xtemp;
    nuc_eos_lookup(xp, xeps, xtemp, xcsq, xrho, xye);

    *vacuum_rho_reset = xrho;
    *vacuum_temp_reset = xtemp;
    *vacuum_ye_reset = 1.2 * eos_yemin;
    *vacuum_P_reset = xp;
    *vacuum_tau_reset = xrho*xeps;

    *vacuum_tau = min(*vacuum_tau, *vacuum_tau_reset);

    *ye_maximum = eos_yemax;


    *minTableEnergy = exp(eos_epsmin) - energy_shift;
    *minTableTemperature = eos_tempmin;

    std::cout<<"vacuum_tau_reset = "<<*vacuum_tau_reset<<std::endl;
    std::cout<<"vacuum_rho_reset = "<<*vacuum_rho_reset<<std::endl;
    std::cout<<"vacuum_P_reset = "<<*vacuum_P_reset<<std::endl;
    std::cout<<"vacuum_temp_reset = "<<*vacuum_temp_reset<<std::endl;
    std::cout<<"vacuum_ye_reset = "<<*vacuum_ye_reset<<std::endl;
    std::cout<<"ye_maximum = "<<*ye_maximum<<std::endl;
    std::cout<<"minTableEnergy = "<<*minTableEnergy<<std::endl;
    std::cout<<"minTableTemperature = "<<*minTableTemperature<<std::endl;
    std::cout<<"****************ExternalInitialData********************"<<std::endl;

    xtemp = eos_tempmin;
    xeps = 0;
    if (*initialYe < 0) {
      xye = ye_beta_equil(xrho, xtemp);
    }
    else {
      xye = *initialYe;
    }

    nuc_eos_lookup(xp, xeps, xtemp, xcsq, xrho, xye);

    double vac_P    = xp;
    double vac_eps  = xeps;
    double vac_ye   = xye;
    double vac_temp = xtemp;
    double vac_csq  = xcsq;

    for (int k=0; k<nz; k++) {
      for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
           int i3D = i + nx*(j  + ny*k);
          xeps = 0;

         // This is the vacuum state. 
          if (lessEq(rho[i3D], eos_rhomin)) {
            rho[i3D]         = eos_rhomin;
            temperature[i3D] = eos_tempmin;
            Ye[i3D]          = vac_ye;
            P[i3D]           = eos_pressmin ;
            eps[i3D]         = eos_epsmin;
            sqcs[i3D]        = vac_csq;
            vx[i3D]          = 0.0;
            vy[i3D]          = 0.0;
            vy[i3D]          = 0.0;  
          }
          else {
            temperature[i3D] = *initialTemperature;
            xtemp = *initialTemperature;
            xrho = rho[i3D];

            if (*initialYe < 0) {
              xye = ye_beta_equil(xrho, xtemp);
            }
            else {
              xye = *initialYe;
            }
            Ye[i3D] = xye;
            nuc_eos_lookup(xp, xeps, xtemp, xcsq, xrho, xye);
            P[i3D]   = xp;
            eps[i3D] = xeps;
            temperature[i3D] = xtemp;
            sqcs[i3D] = xcsq;
            

          }
        }
      }
    }
  }

  magneticField(P, Sx, Sy, Sz, Bx, By, Bz, gxx, gxy, gxz, gyy, gyz, gzz, *tov_angle, *tov_Asize, *tov_facp, *vacuum_rho, init_domain_x, init_domain_y, init_domain_z, init_patch_x, init_patch_y, init_patch_z, d_ghost_width, dx, nx, ny, nz);

  adm2bssn(DYf, Ye, eps, tau, rho, D, vx, vy, vz, Sx, Sy, Sz, P, Bx, By, Bz, gxx, gxy, gxz, gyy, gyz, gzz, Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, Alpha, chi, trK, theta, Gamh_x, Gamh_y, Gamh_z, *Gamma, Betaux, Betauy, Betauz, dx, nx, ny, nz);

}



/*
!----------------------------------------------------------------------
!
!  This routine calculates the sound speed for the nuc eos.
!  A table look up is done.
!
!----------------------------------------------------------------------
*/
void SAMRAI_External_Data::nuc_eos_lookup(double& P, double& eps, double& temp, double& csq, double rho, double ye) {

  int   rc, lrc;
  double    xrho, xtemp, xye, xenr, xprs, xcs2, rfeps;


//       G = 6.67259d-8
//       M = 1.99d33
//       c = 2.99792485d10

  // Set up variables for the table and convert from Geo to CGS
  // units
  //xrho = (G*M/(c*c))**(-3)*M*v(VF_RHO)
  xrho = rho;
  xye  = ye;

  xenr   = eps;
  xtemp  = temp;

  rfeps = 1.0E-10;
  //rfeps = 1.0d-5

  // initialize intent(out) variables
  xprs = 0.0;
  xcs2 = 0.0;


  int keyerr = nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xcs2,rfeps);

  if ( keyerr == 667 ) {
    std::cout <<"nuc_eos_lookup: xrho = "<<xrho<< std::endl;
    std::cout <<"nuc_eos_lookup: xtemp = "<<xtemp<< std::endl;
    std::cout <<"nuc_eos_lookup: xye = "<<xye<< std::endl;
    std::cout <<"nuc_eos_lookup: xenr = "<<xenr<< std::endl;
    keyerr = nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xcs2,rfeps);
    if ( keyerr != 0) {
      std::cout <<"nuc_eos_lookup:  return from nuc_eos_short keyerr = "<<keyerr<< std::endl;
      std::exit(-1);
    }
  }

  //Convert from CGS to Geo units and set up returned quantities

  //P    = xprs/((G*M/(c*c))**(-3)*M*c*c)
  //eps  = xenr/(c*c)
  P    = xprs;
  eps  = xenr;
  temp = xtemp;
  //Divide by h/rho to convert csq as Ott says
  // See Eqn. B.1 in Otts thesis for this.  Note that Eqn B.1 has
  // a typo. The equation in the main body of the dissertation is
  // correct.
  //csq  = xcs2/(1.0d0 + xenr + xprs/xrho)**(2)/(c*c)
  csq  = xcs2/(1.0 + eps + P/rho);
}

int SAMRAI_External_Data::nuc_eos_short(double& xrho, double& xtemp, double xye, double& xenr, double& xprs, double& xcs2, double xrfeps) {

  int keyerr = 0;
  int keyerrt;

  // local variables
  double lr,lt,y,xx,xeps,leps,xs,xpressure;
  double lr0,lt0,y0,leps0,rfeps;
  double d1,d2,d3;
  double lt2, xgam, xcs2_geom; //dwn
  int for_isnan;  //dwn
  double ff[8];
  xenr = 0;

  // slightly modified by DWN to print out rho.
  if (xrho > eos_rhomax) {
    std::cout <<"nuc_eos_short: xrho > eos_rhomax: "<<xrho<<" > "<<eos_rhomax<< std::endl;
    std::exit(-1);
  }

  if(xrho<eos_rhomin*1.2) {
    std::cout <<"nuc_eos_short: xrho < eos_rhomin*1.2: "<<xrho<<" < "<<eos_rhomin*1.2<< std::endl;
    std::exit(-1);
  }

  if(xye > eos_yemax) {
    std::cout <<"nuc_eos_short: xye > eos_yemax: xye, eos_yemax ="<<xye<<" "<<eos_yemax<< std::endl;
  }

  if(xye < eos_yemin) {
     if(xye+1e-5 < eos_yemin) {
      std::cout <<"nuc_eos_short: xye < eos_yemin: xye, eos_yemin ="<<xye<<" "<<eos_yemin<< std::endl;
      std::cout <<"nuc_eos: ye < yemin"<< std::endl;
      std::exit(-1);
     }
  }

  if(xtemp > eos_tempmax) {
    std::cout <<"nuc_eos: temp > tempmax: "<<xtemp<<" > "<<eos_tempmax<< std::endl;
    std::exit(-1);
  }
   
  if(xtemp < eos_tempmin) {
    std::cout <<"nuc_eos: temp < tempmin:"<<xtemp<<" < "<<eos_tempmin<< std::endl;
    std::exit(-1);
  }
  //DWN check on rfeps
  rfeps = max(1.0E-15, xrfeps);

  lr = log(xrho);
  lt = log(xtemp);
  y = xye;
  if(xye > eos_yemax) y = 1.0;
  if(xye < eos_yemin) y = eos_yemin;
  xeps = xenr + energy_shift;
  leps = log(xeps);

  keyerr = 0;


  // have rho,temp,ye; proceed:
  double lrpar[1];
  lrpar[0] = lr;
  double ltpar[1];
  ltpar[0] = lt;
  double ypar[1];
  ypar[0] = y;
  findall_short(lrpar,ltpar,ypar,ff);

  //reset xprs
  xprs = exp(ff[2]);
  if (xprs != xprs) {   // DWN
    std::cout <<"P has a NAN from findall_short. ff[2] = "<<ff[2]<< std::endl;
    std::cout <<lr<<" "<<lt<<" "<<y<<" "<<xrho<< std::endl;
    std::exit(-1);
  }

  //reset xenr
  xenr = exp(ff[1]) - energy_shift;

  xcs2 = ff[0];


  return keyerr;
}

double SAMRAI_External_Data::ye_beta_equil(double rho_geom, double temp) {

  double tol = 1.0E-7;
  double munu;
  double xye;

  double ye_min = yes[0];
  double ye_max = yes[nye-1];

  double rho0 = rho_geom;

  int lrc = be_zbrent(xye, ye_min, ye_max, tol, rho0, temp);
  if (lrc > 0) {
    // solution found
    return xye;
  } else {
    // Find the minimum value of mu_nu as a function of Ye for a given rho0 
    // and T.  Use this value when mu_nu can not be set to zero.
    return find_min_abs_munu(munu, rho0, temp);
  }
}

int SAMRAI_External_Data::be_zbrent(double& x, double x1, double x2, double tol, double rho, double temp) {
  
  int    iter, rc;
  double a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm;
  int itmax = 100;
  double delta = 3.0E-16;

  // Check the endpoints.  The root must be bracketed.
  a = x1;
  b = x2;

  rc = be_eqs(fa, a, rho, temp);
  if (rc < 0) {
    std::cout <<"nuceos_zbrent: failed at endpoint check for a = "<<a<< std::endl;
    return -1;
  }
  if (fabs(fa) < delta) {
    // solution
    x = a;
    return 1;
  }

  rc = be_eqs(fb, b, rho, temp);
  if (rc < 0) {
    std::cout <<"be_zbrent: failed at endpoint check for b = "<<b<< std::endl;
    return -1;
  }
  if (fabs(fb) < delta) {
    // solution
    x = b;
    return 1;
  }

  if(fb*fa > 0.0) {
    return -1;
  }
  fc = fb;
  for (iter = 0; iter < itmax; iter++) {
    if (fb*fc >= 0.0) {
      c  = a;
      fc = fa;
      d  = b-a;
      e  = d;
    }
    if(fabs(fc) < fabs(fb)) {
      a  = b;
      b  = c;
      c  = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0*delta*fabs(b) + 0.5*tol;
    xm = 0.5*(c-b);
    if (fabs(xm) <= tol1 || fabs(fb) < delta) {
      // Success, root found.
      x = b;
      return iter;
    }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb/fa;
      if (a == c) {
        p = 2.*xm*s;
        q = 1. - s;
      } else {
        q = fa/fc;
        r = fb/fc;
        p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if(p > 0.0) q = -q;
      p = fabs(p);
      if (2.0*p < min(3.0*xm*q-fabs(tol1*q),fabs(e*q))) {
        e = d;
        d = p/q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if(fabs(d) > tol1) {
      b = b+d;
    } else {
      double tol1_t = fabs(tol1);
      if (xm < 0) {
        tol1_t = -tol1_t;
      }
      b = b + tol1_t;
    }

    rc = be_eqs(fb, b, rho, temp);
    if (rc < 0) {
      std::cout <<"be_zbrent: call to be_eqs failed in loop for b = "<<b<< std::endl;
    }
  }

  // Solver failed.
  std::cout <<"zbrent exceeding maximum iterations."<< std::endl;

  return -1;
}


double SAMRAI_External_Data::find_min_abs_munu(double& tol, double rho, double temp) {

  double yemin = -1.0;
  double fmin = 1.0E10;
  double f;

  for (int i = 0; i < nye; i++) {
    double xye = yes[i];

    int rc = be_eqs(f, xye, rho, temp);
    if (rc < 0) {
      std::cout <<"find_min_abs_munu: failed"<< std::endl;
      std::exit(-1);
    }
    if (fabs(f) < fmin) {
      fmin = fabs(f);
      yemin = xye;
    }
  }

  tol = fmin;

  return yemin;
}

double SAMRAI_External_Data::be_eqs(double& f, double x, double rho, double temp) {


  double xrho, xtemp, xye,xmu_e,xmu_n,xmu_p;

  xrho = rho;
  xtemp = temp;
  xye = x;

  check_eostable_input(xrho, xtemp, xye);
  nuc_eos_mu(xrho,xtemp,xye,xmu_e,xmu_n,xmu_p);

  f =  xmu_p + xmu_e - xmu_n;

  return 1;
}

void SAMRAI_External_Data::check_eostable_input(double xrho,double xtemp,double& xye) {
  int e = 0;

  if (xrho < eos_rhomin) {
    std::cout <<"check_eostable_input: xrho < eos_rhomin: "<<xrho<<" "<<eos_rhomin<< std::endl;
    std::exit(-1);
  }

  // only check the temperature limits if keytemp != 0.
  if (xtemp < eos_tempmin) {
    std::cout <<"check_eostable_input: xtemp < eos_tempmin: "<<xtemp<<" "<<eos_tempmin<< std::endl;
    std::exit(-1);
  }

  if (xtemp > eos_tempmax) {
    std::cout <<"check_eostable_input: xtemp > eos_tempmax: "<<xtemp<<" "<<eos_tempmax<< std::endl;
    std::exit(-1);
  }

  if (xye < eos_yemin) {
   if (xye - eos_yemin < -1.0e-8) {
    std::cout <<"check_eostable_input: xye < eos_yemin: "<<xye<<" "<<eos_yemin<< std::endl;
    std::exit(-1);
   }
   xye = eos_yemin;
  }

  if (xye > eos_yemax) {
    std::cout <<"check_eostable_input: xye > eos_yemax: "<<xye<<" "<<eos_yemax<< std::endl;
    std::exit(-1);
  }

}

void SAMRAI_External_Data::nuc_eos_mu(double& xrho,double& xtemp,double xye,double& xmu_e,double& xmu_n,double& xmu_p) {


  double lr,lt,y;
  double ff[3];
  
  if (xrho > eos_rhomax) {
    std::cout <<"nuc_eos_full: rho > rhomax: rho, rhomax = "<<xrho<<" "<<eos_rhomax<< std::endl;
    std::exit(-1);
  }

  if(xrho < eos_rhomin) {
    std::cout <<"nuc_eos: rho < rhomin"<< std::endl;
    std::exit(-1);    
  }

  if(xye > eos_yemax) {
    std::cout <<"nuc_eos_full: xye > eos_yemax: xye, eos_yemax ="<<xye<<" "<<eos_yemax<< std::endl;
  }

  if(xye < eos_yemin) {
     if(xye+1e-5 < eos_yemin) {
      std::cout <<"nuc_eos_full: xye < eos_yemin: xye, eos_yemin ="<<xye<<" "<<eos_yemin<< std::endl;
      std::cout <<"nuc_eos: ye < yemin"<< std::endl;
      std::exit(-1);
     }
  }

  if(xtemp > eos_tempmax) {
    std::cout <<"nuc_eos: temp > tempmax"<< std::endl;
    std::exit(-1);
  }
   
  if(xtemp < eos_tempmin) {
    std::cout <<"nuc_eos: temp < tempmin"<< std::endl;
    std::exit(-1);
  }

  lr = log(xrho);
  lt = log(xtemp);
  y = xye;
  if(xye > eos_yemax) y = 1.0;
  if(xye < eos_yemin) y = eos_yemin;

  // have rho,T,ye; proceed:
  double lrpar[1];
  lrpar[0] = lr;
  double ltpar[1];
  ltpar[0] = lt;
  double ypar[1];
  ypar[0] = y;
  findmus(lrpar,ltpar,ypar,ff);

  
// chemical potentials
  xmu_e = ff[0];
  xmu_p = ff[1];
  xmu_n = ff[2];

}


void SAMRAI_External_Data::findmus(double* lr,double* lt, double* y, double* ff) {
  double d1,d2,d3;
  // Ewald's interpolator           
  int ivs[3];
  for (int i = 0; i < 3; i++) {
    ivs[i] = i + 3;
  }
  intp3d(lr,lt,y,ff,1,alltables,ivs, 3, nrho,ntemp,nye,logrho,logtemp,yes, d1, d2, d3);
}

void SAMRAI_External_Data::findall_short(double* lr,double* lt, double* y, double* ff) {
  double d1,d2,d3;
  // Ewald's interpolator           
  int ivs[3];
  for (int i = 0; i < 3; i++) {
    ivs[i] = i;
  }
  intp3d(lr,lt,y,ff,1,alltables,ivs, 3, nrho,ntemp,nye,logrho,logtemp,yes, d1, d2, d3);
}

/*
c                                                          
c---------------------------------------------------------------------
c
c     purpose: interpolation of a function of three variables in an
c              equidistant(!!!) table.
c
c     method:  8-point Lagrange linear interpolation formula          
c
c     x        input vector of first  variable
c     y        input vector of second variable
c     z        input vector of third  variable
c
c     f        output vector of interpolated function values
c
c     kt       vector length of input and output vectors
c
c     ft       3d array of tabulated function values
c     nx       x-dimension of table
c     ny       y-dimension of table
c     nz       z-dimension of table
c     xt       vector of x-coordinates of table
c     yt       vector of y-coordinates of table
c     zt       vector of z-coordinates of table
c
c     d1       centered derivative of ft with respect to x
c     d2       centered derivative of ft with respect to y
c     d3       centered derivative of ft with respect to z
c     Note that d? only make sense when intp3d is called with kt=1
c---------------------------------------------------------------------
c
c 
*/

void SAMRAI_External_Data::intp3d(double* x, double* y, double* z, double* f, int kt, double* ft, int* ivs, int nivs, int nx, int ny, int nz, double* xt, double* yt, double* zt, double& d1, double& d2, double& d3) {

  //------  determine spacing parameters of (equidistant!!!) table

  double dx    = (xt[nx - 1] - xt[0]) / (nx-1);
  double dy    = (yt[ny - 1] - yt[0]) / (ny-1);
  double dz    = (zt[nz - 1] - zt[0]) / (nz-1);

  double dxi   = 1. / dx;
  double dyi   = 1. / dy;
  double dzi   = 1. / dz;

  double dxyi  = dxi * dyi;
  double dxzi  = dxi * dzi;
  double dyzi  = dyi * dzi;

  double dxyzi = dxi * dyi * dzi;


  double fh[8];
  int idx[8];

  //------- loop over all points to be interpolated

  for (int n = 0; n < kt; n++) {

  //------- determine location in (equidistant!!!) table 
                                                                  
    int ix = 1 + int( (x[n] - xt[0] - 1.e-10) * dxi );
    int iy = 1 + int( (y[n] - yt[0] - 1.e-10) * dyi );
    int iz = 1 + int( (z[n] - zt[0] - 1.e-10) * dzi );
                                                     
    ix = MAX( 1, MIN( ix, nx - 1 ) );
    iy = MAX( 1, MIN( iy, ny - 1 ) );
    iz = MAX( 1, MIN( iz, nz - 1 ) );


  //------- set-up auxiliary arrays for Lagrange interpolation

    double delx = xt[ix] - x[n];
    double dely = yt[iy] - y[n];
    double delz = zt[iz] - z[n];

    idx[0] = NTABLES*(ix + nx*(iy + ny*iz));
    idx[1] = NTABLES*((ix-1) + nx*(iy + ny*iz));
    idx[2] = NTABLES*(ix + nx*((iy-1) + ny*iz));
    idx[3] = NTABLES*(ix + nx*(iy + ny*(iz-1)));
    idx[4] = NTABLES*((ix-1) + nx*((iy-1) + ny*iz));
    idx[5] = NTABLES*((ix-1) + nx*(iy + ny*(iz-1)));
    idx[6] = NTABLES*(ix + nx*((iy-1) + ny*(iz-1)));
    idx[7] = NTABLES*((ix-1) + nx*((iy-1) + ny*(iz-1)));



    for (int i = 0; i < nivs; i++) {
      fh[0] = ft[ivs[i]+idx[0]];
      fh[1] = ft[ivs[i]+idx[1]];
      fh[2] = ft[ivs[i]+idx[2]];
      fh[3] = ft[ivs[i]+idx[3]];
      fh[4] = ft[ivs[i]+idx[4]];
      fh[5] = ft[ivs[i]+idx[5]];
      fh[6] = ft[ivs[i]+idx[6]];
      fh[7] = ft[ivs[i]+idx[7]];


      //------ set up coefficients of the interpolation polynomial and 
      //       evaluate function values 
                                                      
      double a1 = fh[0];
      double a2 = dxi   * ( fh[1] - fh[0] );
      double a3 = dyi   * ( fh[2] - fh[0] );
      double a4 = dzi   * ( fh[3] - fh[0] );
      double a5 = dxyi  * ( fh[4] - fh[1] - fh[2] + fh[0] );
      double a6 = dxzi  * ( fh[5] - fh[1] - fh[3] + fh[0] );
      double a7 = dyzi  * ( fh[6] - fh[2] - fh[3] + fh[0] );
      double a8 = dxyzi * ( fh[7] - fh[0] + fh[1] + fh[2] + fh[3] - fh[4] - fh[5] - fh[6] );

      d1 = -a2;
      d2 = -a3;
      d3 = -a4;

      f[i + n * nivs]  = a1 +  a2 * delx +  a3 * dely +  a4 * delz +  a5 * delx * dely +  a6 * delx * delz +  a7 * dely * delz +  a8 * delx * dely * delz;
    }

  }                                                               
}                                                          

void SAMRAI_External_Data::nuc_eos_C_ReadTable(char* nuceos_table_name)
{
//  double temp0, temp1;
//  double energy_shift;
//
//  double eos_rhomax, eos_rhomin;
//  double eos_tempmin, eos_tempmax;
//  double eos_yemin, eos_yemax;
//  double eos_epsmin, eos_epsmax;
//  
//  double c2p_tempmin;
//  double c2p_tempmax;
//
//  int nrho;
//  int ntemp;
//  int nye;
//
//  double * restrict alltables;
//  double * restrict epstable;
//  double * restrict logrho;
//  double * restrict logtemp;
//  double dlintemp, dlintempi;
//  double drholintempi;
//  double dlintempyei;
//  double drholintempyei;
//  double * restrict yes;
//  double dtemp, dtempi;
//  double drho, drhoi;
//  double dye, dyei;
//  double drhotempi;
//  double drhoyei;
//  double dtempyei;
//  double drhotempyei;

  printf("*******************************\n");
  printf("Reading nuc_eos table file:\n");
  printf("%s\n",nuceos_table_name);
  printf("*******************************\n");

  hid_t file;
  file = H5Fopen(nuceos_table_name, H5F_ACC_RDONLY, H5P_DEFAULT);

// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_EOS_HDF5(NAME,VAR,TYPE,MEM)                                      \
  do {                                                                        \
    hid_t dataset;                                                            \
    dataset = H5Dopen(file, NAME, H5P_DEFAULT);                                \
    H5Dread(dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR);       \
    H5Dclose(dataset);                                            \
  } while (0)
// The second reads a given variable into a hyperslab of the alltables_temp array
#define READ_EOSTABLE_HDF5(NAME,OFF)                                     \
  do {                                                                   \
    hsize_t offset[2]     = {OFF,0};                                     \
    H5Sselect_hyperslab(mem3, H5S_SELECT_SET, offset, NULL, var3, NULL); \
    READ_EOS_HDF5(NAME,alltables_temp,H5T_NATIVE_DOUBLE,mem3);           \
  } while (0)


  // Read logrho, logtemp and ye
  int rank;
  hid_t dataspace;
  std::string coordName;
  hid_t dataset;
  hsize_t* dims;
  dataset = H5Dopen(file, "coord0", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  rank = H5Sget_simple_extent_ndims(dataspace);
  dims = new hsize_t[rank];
  H5Sget_simple_extent_dims(dataspace, dims, dims);
  logrho = new double[dims[0]];
  nrho = dims[0];
  READ_EOS_HDF5("coord0", logrho, H5T_NATIVE_DOUBLE, H5S_ALL);
  dataset = H5Dopen(file, "coord1", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  rank = H5Sget_simple_extent_ndims(dataspace);
  dims = new hsize_t[rank];
  H5Sget_simple_extent_dims(dataspace, dims, dims);
  logtemp = new double[dims[0]];
  ntemp = dims[0];
  READ_EOS_HDF5("coord1", logtemp, H5T_NATIVE_DOUBLE, H5S_ALL);
  dataset = H5Dopen(file, "coord2", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);
  rank = H5Sget_simple_extent_ndims(dataspace);
  dims = new hsize_t[rank];
  H5Sget_simple_extent_dims(dataspace, dims, dims);
  yes = new double[dims[0]];
  nye = dims[0];
  READ_EOS_HDF5("coord2", yes, H5T_NATIVE_DOUBLE, H5S_ALL);

  // Allocate memory for tables
  double* alltables_temp;
  if (!(alltables_temp = (double*)malloc(nrho * ntemp * nye * NTABLES * sizeof(double)))) {
    printf("Cannot allocate memory for EOS table");
  }

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[2] = {NTABLES, (hsize_t)nrho * ntemp * nye};
  hsize_t var3[2]       = { 1, (hsize_t)nrho * ntemp * nye};
  hid_t mem3 =  H5Screate_simple(2, table_dims, NULL);

  // Read alltables_temp
  READ_EOSTABLE_HDF5("var0",  0);
  READ_EOSTABLE_HDF5("var1",  1);
  READ_EOSTABLE_HDF5("var2",  2);
  READ_EOSTABLE_HDF5("var3",  3);
  READ_EOSTABLE_HDF5("var4",  4);
  READ_EOSTABLE_HDF5("var5",  5);


  H5Sclose(mem3);
  H5Fclose(file);

  // change ordering of alltables array so that
  // the table kind is the fastest changing index
  if (!(alltables = (double*)malloc(nrho * ntemp * nye * NTABLES 
            * sizeof(double)))) {
    printf("Cannot allocate memory for EOS table");
  }
  int idxEPS = 1;
  eos_epsmin = alltables_temp[nrho*(ntemp*(nye*idxEPS))];
  for(int iv = 0;iv<NTABLES;iv++) 
    for(int k = 0; k<nye;k++) 
      for(int j = 0; j<ntemp; j++) 
  for(int i = 0; i<nrho; i++) {
    int indold = i + nrho*(j + ntemp*(k + nye*iv));
    int indnew = iv + NTABLES*(i + nrho*(j + ntemp*k));
    alltables[indnew] = alltables_temp[indold];
    //Get epsmin
    if (iv == idxEPS && alltables[indnew] < eos_epsmin) {
      eos_epsmin = alltables[indnew];
    }
  }

  // free memory of temporary array
  free(alltables_temp);

  
  eos_rhomin = exp(logrho[0]);
  eos_rhomax = exp(logrho[nrho-1]);

  eos_yemax = yes[nye-1];
  eos_yemin = yes[0];

  eos_tempmin = exp(logtemp[0]);
  eos_tempmax = exp(logtemp[ntemp-1]);


}


void SAMRAI_External_Data::magneticField(double *P,   double *Sx,     double *Sy,     double *Sz, double *Bx,     double *By,     double *Bz, 
              double *gxx,   double *gxy,    double *gxz,
              double *gyy,   double *gyz,    double *gzz,
              double tov_angle, double tov_Asize, double tov_facp, double vacuum, double init_domain_x, double init_domain_y, double init_domain_z, int init_patch_x, int init_patch_y, int init_patch_z, int d_ghost_width, const double* dxA, int nx,        int ny,         int nz) {

  double dx = dxA[0];
  double dy = dxA[1];
  double dz = dxA[2];

  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        int i3D = i + nx*(j  + ny*k);

        // convert position from cactus units to Lorene units (10Km)
        double x_fd = (init_domain_x + ((i + init_patch_x) - d_ghost_width) * dxA[0])*CUto10KM;
        double y_fd = (init_domain_y + ((j + init_patch_y) - d_ghost_width) * dxA[1])*CUto10KM;
        double z_fd = (init_domain_z + ((k + init_patch_z) - d_ghost_width) * dxA[2])*CUto10KM;

        double myz     = -y_fd*sin(tov_angle) + z_fd*cos(tov_angle);
        double myy     =  y_fd*cos(tov_angle) + z_fd*sin(tov_angle);
        double mytmp   =  tov_Asize * max( P[i3D]-tov_facp*vacuum, 0.0);
        Sx[i3D]  =                 ( -myy  * mytmp );
        Sy[i3D]  =  cos(tov_angle) * (  x_fd  * mytmp );
        Sz[i3D]  =  sin(tov_angle) * (  x_fd  * mytmp );
      }
    }
  }

  for (int k=3; k<nz-3; k++) {
    for (int j=3; j<ny-3; j++) {
      for (int i=3; i<nx-3; i++) {
        int i3D = i + nx*(j  + ny*k);

        int i3Dip1 = i + 1 + nx*(j  + ny*k);
        int i3Dip2 = i + 2 + nx*(j  + ny*k);
        int i3Dip3 = i + 3 + nx*(j  + ny*k);
        int i3Dim1 = i - 1 + nx*(j  + ny*k);
        int i3Dim2 = i - 2 + nx*(j  + ny*k);
        int i3Dim3 = i - 3 + nx*(j  + ny*k);

        int i3Djp1 = i + nx*(j + 1 + ny*k);
        int i3Djp2 = i + nx*(j + 2 + ny*k);
        int i3Djp3 = i + nx*(j + 3 + ny*k);
        int i3Djm1 = i + nx*(j - 1 + ny*k);
        int i3Djm2 = i + nx*(j - 2 + ny*k);
        int i3Djm3 = i + nx*(j - 3 + ny*k);

        int i3Dkp1 = i + nx*(j  + ny*(k + 1));
        int i3Dkp2 = i + nx*(j  + ny*(k + 2));
        int i3Dkp3 = i + nx*(j  + ny*(k + 3));
        int i3Dkm1 = i + nx*(j  + ny*(k - 1));
        int i3Dkm2 = i + nx*(j  + ny*(k - 2));
        int i3Dkm3 = i + nx*(j  + ny*(k - 3));


        double dxAx = CD6(Sx, i3Dip3, i3Dip2, i3Dip1, i3Dim1, i3Dim2, i3Dim3, dx);
        double dyAx = CD6(Sx, i3Djp3, i3Djp2, i3Djp1, i3Djm1, i3Djm2, i3Djm3, dy);
        double dzAx = CD6(Sx, i3Dkp3, i3Dkp2, i3Dkp1, i3Dkm1, i3Dkm2, i3Dkm3, dz);

        double dxAy = CD6(Sy, i3Dip3, i3Dip2, i3Dip1, i3Dim1, i3Dim2, i3Dim3, dx);
        double dyAy = CD6(Sy, i3Djp3, i3Djp2, i3Djp1, i3Djm1, i3Djm2, i3Djm3, dy);
        double dzAy = CD6(Sy, i3Dkp3, i3Dkp2, i3Dkp1, i3Dkm1, i3Dkm2, i3Dkm3, dz);

        double dxAz = CD6(Sz, i3Dip3, i3Dip2, i3Dip1, i3Dim1, i3Dim2, i3Dim3, dx);
        double dyAz = CD6(Sz, i3Djp3, i3Djp2, i3Djp1, i3Djm1, i3Djm2, i3Djm3, dy);
        double dzAz = CD6(Sz, i3Dkp3, i3Dkp2, i3Dkp1, i3Dkm1, i3Dkm2, i3Dkm3, dz);        

        double gtd_xx = gxx[i3D];
        double gtd_xy = gxy[i3D];
        double gtd_xz = gxz[i3D];
        double gtd_yy = gyy[i3D];
        double gtd_yz = gyz[i3D];
        double gtd_zz = gzz[i3D];
        
        double sdetg_pt = sqrt(gtd_xx*gtd_yy*gtd_zz-gtd_xx*gtd_yz*gtd_yz-gtd_xy*gtd_xy*gtd_zz+2*gtd_xy*gtd_xz*gtd_yz-gtd_xz*gtd_xz*gtd_yy);

        Bx[i3D]   = (-1./sdetg_pt) *(-dzAy + dyAz);
        By[i3D]   = (-1./sdetg_pt) *( dzAx - dxAz);
        Bz[i3D]   = (-1./sdetg_pt) *( dxAy - dyAx);
      }
    }
  }

}


void SAMRAI_External_Data::adm2bssn(double *DYf, double *Ye, double *eps, double *tau, double *rho, double *D,   double *vx,     double *vy,     double *vz,   double *Sx,     double *Sy,     double *Sz,
              double *P,     double *Bx,     double *By,     double *Bz,
              double *gxx,   double *gxy,    double *gxz,
              double *gyy,   double *gyz,    double *gzz,
              double *Kxx,   double *Kxy,    double *Kxz,
              double *Kyy,   double *Kyz,    double *Kzz,
              double *Alpha, double *chi, double *trK, double *theta,
              double *Gamh_x,     double *Gamh_y,     double *Gamh_z, double Gamma, 
              double *Betaux,double *Betauy, double *Betauz, const double* dxA,
              int nx,        int ny,         int nz) {

  int nd = nx*ny*nz;
  double dx = dxA[0];
  double dy = dxA[1];
  double dz = dxA[2];
      
  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        int i3D = i + nx*(j  + ny*k);
        double gd_xx = gxx[i3D];
        double gd_xy = gxy[i3D];
        double gd_xz = gxz[i3D];
        double gd_yy = gyy[i3D];
        double gd_yz = gyz[i3D];
        double gd_zz = gzz[i3D];
    
        double vux = vx[i3D];
        double vuy = vy[i3D];
        double vuz = vz[i3D];

        double Kd_xx = Kxx[i3D];
        double Kd_xy = Kxy[i3D];
        double Kd_xz = Kxz[i3D];
        double Kd_yy = Kyy[i3D];
        double Kd_yz = Kyz[i3D];
        double Kd_zz = Kzz[i3D];
   
        double t1 = gd_xx;
        double t2 = gd_yy;
        double t4 = gd_zz;
        double t6 = gd_yz;
        double t7 = t6*t6;
        double t9 = gd_xy;
        double t10 = t9*t9;
        double t12 = gd_xz;
        double t16 = t12*t12;
        double detgd = t1*t2*t4-t1*t7-t10*t4+2*t9*t12*t6-t16*t2;
        double idetgd = 1.0/detgd;
        double gu_xx = idetgd*(gd_yy*gd_zz-gd_yz*gd_yz);
        double gu_xy = idetgd*(-gd_xy*gd_zz+gd_xz*gd_yz);
        double gu_xz = idetgd*(gd_xy*gd_yz-gd_xz*gd_yy);
        double gu_yy = idetgd*(gd_xx*gd_zz-gd_xz*gd_xz);
        double gu_yz = idetgd*(-gd_xx*gd_yz+gd_xy*gd_xz);
        double gu_zz = idetgd*(gd_xx*gd_yy-gd_xy*gd_xy);

        double chitmp = pow(idetgd, (1.0/3.0));
        chi[i3D] = chitmp;

        gxx[i3D] = chitmp*gd_xx;
        gxy[i3D] = chitmp*gd_xy;
        gxz[i3D] = chitmp*gd_xz;
        gyy[i3D] = chitmp*gd_yy;
        gyz[i3D] = chitmp*gd_yz;
        gzz[i3D] = chitmp*gd_zz;


        double gtd_xx = gxx[i3D];
        double gtd_xy = gxy[i3D];
        double gtd_xz = gxz[i3D];
        double gtd_yy = gyy[i3D];
        double gtd_yz = gyz[i3D];
        double gtd_zz = gzz[i3D];

        double TrK = gu_xx*Kd_xx+gu_yy*Kd_yy+gu_zz*Kd_zz+2*(gu_xy*Kd_xy+gu_xz*Kd_xz+gu_yz*Kd_yz);
        trK[i3D] = TrK;

        Kxx[i3D] = chitmp*(Kd_xx-(1.0/3.0)*gd_xx*TrK);
        Kxy[i3D] = chitmp*(Kd_xy-(1.0/3.0)*gd_xy*TrK);
        Kxz[i3D] = chitmp*(Kd_xz-(1.0/3.0)*gd_xz*TrK);
        Kyy[i3D] = chitmp*(Kd_yy-(1.0/3.0)*gd_yy*TrK);
        Kyz[i3D] = chitmp*(Kd_yz-(1.0/3.0)*gd_yz*TrK);
        Kzz[i3D] = chitmp*(Kd_zz-(1.0/3.0)*gd_zz*TrK);

        //Primitives to conservatives
        vx[i3D] = gd_xx * vux + gd_xy * vuy + gd_xz * vuz;
        vy[i3D] = gd_xy * vux + gd_yy * vuy + gd_yz * vuz;
        vz[i3D] = gd_xz * vux + gd_yz * vuy + gd_zz * vuz;
        double vfd_x = vx[i3D];
        double vfd_y = vy[i3D];
        double vfd_z = vz[i3D];
        double chi_max = MAX(0.0001, chitmp);
        double h = rho[i3D] * (1.0 + eps[i3D]) + P[i3D];
        double sdetg = pow(chi[i3D], (-1.500000000000000));
        double Bu_x = Bx[i3D];
        double Bu_y = By[i3D];
        double Bu_z = Bz[i3D];
        double inv_chi = 1.0 / fabs(chi_max);
        double Bfd_x = inv_chi * gtd_xx * Bx[i3D] + inv_chi * gtd_xy * By[i3D] + inv_chi * gtd_xz * Bz[i3D];
        double Bfd_y = inv_chi * gtd_xy * Bx[i3D] + inv_chi * gtd_yy * By[i3D] + inv_chi * gtd_yz * Bz[i3D];
        double Bfd_z = inv_chi * gtd_xz * Bx[i3D] + inv_chi * gtd_yz * By[i3D] + inv_chi * gtd_zz * Bz[i3D];

        double detgtd = gtd_xx*gtd_yy*gtd_zz-gtd_xx*gtd_yz*gtd_yz-gtd_xy*gtd_xy*gtd_zz+2*gtd_xy*gtd_xz*gtd_yz-gtd_xz*gtd_xz*gtd_yy;
        double idetgtd = 1.0/detgtd;
        double gtu_xx = idetgtd*(gtd_yy*gtd_zz-gtd_yz*gtd_yz);
        double gtu_xy = idetgtd*(-gtd_xy*gtd_zz+gtd_xz*gtd_yz);
        double gtu_xz = idetgtd*(gtd_xy*gtd_yz-gtd_xz*gtd_yy);
        double gtu_yy = idetgtd*(gtd_xx*gtd_zz-gtd_xz*gtd_xz);
        double gtu_yz = idetgtd*(-gtd_xx*gtd_yz+gtd_xy*gtd_xz);
        double gtu_zz = idetgtd*(gtd_xx*gtd_yy-gtd_xy*gtd_xy);

        double W = 1.0 / sqrt(1.0 - (vfd_x * vux + vfd_y * vuy + vfd_z * vuz));
        double Bd_x = Bfd_x;
        double Bd_y = Bfd_y;
        double Bd_z = Bfd_z;
        D[i3D] = sdetg * rho[i3D] * W;
        DYf[i3D] = D[i3D] * Ye[i3D];

        Sx[i3D] = sdetg * ((h * W * W + (Bd_x * Bu_x + Bd_y * Bu_y + Bd_z * Bu_z)) * vfd_x - (Bu_x * vfd_x + Bu_y * vfd_y + Bu_z * vfd_z) * Bd_x);
        Sy[i3D] = sdetg * ((h * W * W + (Bd_x * Bu_x + Bd_y * Bu_y + Bd_z * Bu_z)) * vfd_y - (Bu_x * vfd_x + Bu_y * vfd_y + Bu_z * vfd_z) * Bd_y);
        Sz[i3D] = sdetg * ((h * W * W + (Bd_x * Bu_x + Bd_y * Bu_y + Bd_z * Bu_z)) * vfd_z - (Bu_x * vfd_x + Bu_y * vfd_y + Bu_z * vfd_z) * Bd_z);
        tau[i3D] = sdetg * (h * W * W + (-P[i3D]) + (-rho[i3D] * W) + (Bd_x * Bu_x + Bd_y * Bu_y + Bd_z * Bu_z) + (-((Bd_x * Bu_x + Bd_y * Bu_y + Bd_z * Bu_z) / (W * W) + (Bu_x * vfd_x + Bu_y * vfd_y + Bu_z * vfd_z) * (Bu_x * vfd_x + Bu_y * vfd_y + Bu_z * vfd_z)) / 2.0));


        Bx[i3D] = sdetg * Bx[i3D];
        By[i3D] = sdetg * By[i3D];
        Bz[i3D] = sdetg * Bz[i3D];
      }
    }
  }

  double *dx_gtd11 = new double[nd];
  double *dx_gtd12 = new double[nd];
  double *dx_gtd13 = new double[nd];
  double *dx_gtd22 = new double[nd];
  double *dx_gtd23 = new double[nd];
  double *dx_gtd33 = new double[nd];

  double *dy_gtd11 = new double[nd];
  double *dy_gtd12 = new double[nd];
  double *dy_gtd13 = new double[nd];
  double *dy_gtd22 = new double[nd];
  double *dy_gtd23 = new double[nd];
  double *dy_gtd33 = new double[nd];

  double *dz_gtd11 = new double[nd];
  double *dz_gtd12 = new double[nd];
  double *dz_gtd13 = new double[nd];
  double *dz_gtd22 = new double[nd];
  double *dz_gtd23 = new double[nd];
  double *dz_gtd33 = new double[nd];

  for (int k=3; k<nz-3; k++) {
    for (int j=3; j<ny-3; j++) {
      for (int i=3; i<nx-3; i++) {
        int i3D = i + nx*(j  + ny*k);


        int i3Dip1 = i + 1 + nx*(j  + ny*k);
        int i3Dip2 = i + 2 + nx*(j  + ny*k);
        int i3Dip3 = i + 3 + nx*(j  + ny*k);
        int i3Dim1 = i - 1 + nx*(j  + ny*k);
        int i3Dim2 = i - 2 + nx*(j  + ny*k);
        int i3Dim3 = i - 3 + nx*(j  + ny*k);

        int i3Djp1 = i + nx*(j + 1 + ny*k);
        int i3Djp2 = i + nx*(j + 2 + ny*k);
        int i3Djp3 = i + nx*(j + 3 + ny*k);
        int i3Djm1 = i + nx*(j - 1 + ny*k);
        int i3Djm2 = i + nx*(j - 2 + ny*k);
        int i3Djm3 = i + nx*(j - 3 + ny*k);

        int i3Dkp1 = i + nx*(j  + ny*(k + 1));
        int i3Dkp2 = i + nx*(j  + ny*(k + 2));
        int i3Dkp3 = i + nx*(j  + ny*(k + 3));
        int i3Dkm1 = i + nx*(j  + ny*(k - 1));
        int i3Dkm2 = i + nx*(j  + ny*(k - 2));
        int i3Dkm3 = i + nx*(j  + ny*(k - 3));


        dx_gtd11[i3D] = CD6(gxx, i3Dip3, i3Dip2, i3Dip1, i3Dim1, i3Dim2, i3Dim3, dx);
        dy_gtd11[i3D] = CD6(gxx, i3Djp3, i3Djp2, i3Djp1, i3Djm1, i3Djm2, i3Djm3, dy);
        dz_gtd11[i3D] = CD6(gxx, i3Dkp3, i3Dkp2, i3Dkp1, i3Dkm1, i3Dkm2, i3Dkm3, dz);

        dx_gtd12[i3D] = CD6(gxy, i3Dip3, i3Dip2, i3Dip1, i3Dim1, i3Dim2, i3Dim3, dx);
        dy_gtd12[i3D] = CD6(gxy, i3Djp3, i3Djp2, i3Djp1, i3Djm1, i3Djm2, i3Djm3, dy);
        dz_gtd12[i3D] = CD6(gxy, i3Dkp3, i3Dkp2, i3Dkp1, i3Dkm1, i3Dkm2, i3Dkm3, dz);

        dx_gtd13[i3D] = CD6(gxz, i3Dip3, i3Dip2, i3Dip1, i3Dim1, i3Dim2, i3Dim3, dx);
        dy_gtd13[i3D] = CD6(gxz, i3Djp3, i3Djp2, i3Djp1, i3Djm1, i3Djm2, i3Djm3, dy);
        dz_gtd13[i3D] = CD6(gxz, i3Dkp3, i3Dkp2, i3Dkp1, i3Dkm1, i3Dkm2, i3Dkm3, dz);

        dx_gtd22[i3D] = CD6(gyy, i3Dip3, i3Dip2, i3Dip1, i3Dim1, i3Dim2, i3Dim3, dx);
        dy_gtd22[i3D] = CD6(gyy, i3Djp3, i3Djp2, i3Djp1, i3Djm1, i3Djm2, i3Djm3, dy);
        dz_gtd22[i3D] = CD6(gyy, i3Dkp3, i3Dkp2, i3Dkp1, i3Dkm1, i3Dkm2, i3Dkm3, dz);

        dx_gtd23[i3D] = CD6(gyz, i3Dip3, i3Dip2, i3Dip1, i3Dim1, i3Dim2, i3Dim3, dx);
        dy_gtd23[i3D] = CD6(gyz, i3Djp3, i3Djp2, i3Djp1, i3Djm1, i3Djm2, i3Djm3, dy);
        dz_gtd23[i3D] = CD6(gyz, i3Dkp3, i3Dkp2, i3Dkp1, i3Dkm1, i3Dkm2, i3Dkm3, dz);

        dx_gtd33[i3D] = CD6(gzz, i3Dip3, i3Dip2, i3Dip1, i3Dim1, i3Dim2, i3Dim3, dx);
        dy_gtd33[i3D] = CD6(gzz, i3Djp3, i3Djp2, i3Djp1, i3Djm1, i3Djm2, i3Djm3, dy);
        dz_gtd33[i3D] = CD6(gzz, i3Dkp3, i3Dkp2, i3Dkp1, i3Dkm1, i3Dkm2, i3Dkm3, dz);
      }
    }
  }

  for (int k=3; k<nz-3; k++) {
    for (int j=3; j<ny-3; j++) {
      for (int i=3; i<nx-3; i++) {
        int i3D = i + nx*(j  + ny*k);
        double gtd_xx = gxx[i3D];
        double gtd_xy = gxy[i3D];
        double gtd_xz = gxz[i3D];
        double gtd_yy = gyy[i3D];
        double gtd_yz = gyz[i3D];
        double gtd_zz = gzz[i3D];
            
        double Dgtd_xxx = dx_gtd11[i3D];
        double Dgtd_xxy = dx_gtd12[i3D];
        double Dgtd_xxz = dx_gtd13[i3D];
        double Dgtd_xyx = Dgtd_xxy;
        double Dgtd_xyy = dx_gtd22[i3D];
        double Dgtd_xyz = dx_gtd23[i3D];
        double Dgtd_xzx = Dgtd_xxz;
        double Dgtd_xzy = Dgtd_xyz;
        double Dgtd_xzz = dx_gtd33[i3D];

        double Dgtd_yxx = dy_gtd11[i3D];
        double Dgtd_yxy = dy_gtd12[i3D];
        double Dgtd_yxz = dy_gtd13[i3D];
        double Dgtd_yyx = Dgtd_yxy;
        double Dgtd_yyy = dy_gtd22[i3D];
        double Dgtd_yyz = dy_gtd23[i3D];
        double Dgtd_yzx = Dgtd_yxz;
        double Dgtd_yzy = Dgtd_yyz;
        double Dgtd_yzz = dy_gtd33[i3D];
    
        double Dgtd_zxx = dz_gtd11[i3D];
        double Dgtd_zxy = dz_gtd12[i3D];
        double Dgtd_zxz = dz_gtd13[i3D];
        double Dgtd_zyx = Dgtd_zxy;
        double Dgtd_zyy = dz_gtd22[i3D];
        double Dgtd_zyz = dz_gtd23[i3D];
        double Dgtd_zzx = Dgtd_zxz;
        double Dgtd_zzy = Dgtd_zyz;
        double Dgtd_zzz = dz_gtd33[i3D];

        double detgtd = gtd_xx*gtd_yy*gtd_zz-gtd_xx*gtd_yz*gtd_yz-gtd_xy*gtd_xy*gtd_zz+2*gtd_xy*gtd_xz*gtd_yz-gtd_xz*gtd_xz*gtd_yy;
        double idetgtd = 1.0/detgtd;
        double gtu_xx = idetgtd*(gtd_yy*gtd_zz-gtd_yz*gtd_yz);
        double gtu_xy = idetgtd*(-gtd_xy*gtd_zz+gtd_xz*gtd_yz);
        double gtu_xz = idetgtd*(gtd_xy*gtd_yz-gtd_xz*gtd_yy);
        double gtu_yy = idetgtd*(gtd_xx*gtd_zz-gtd_xz*gtd_xz);
        double gtu_yz = idetgtd*(-gtd_xx*gtd_yz+gtd_xy*gtd_xz);
        double gtu_zz = idetgtd*(gtd_xx*gtd_yy-gtd_xy*gtd_xy);
        double t3 = Dgtd_yxz+Dgtd_zxy-Dgtd_xyz;
        double t8 = Dgtd_xyz+Dgtd_zxy-Dgtd_yxz;
        double t13 = Dgtd_xyz+Dgtd_yxz-Dgtd_zxy;
        double Ctd_xxx = Dgtd_xxx;
        double Ctd_xxy = Dgtd_yxx;
        double Ctd_xxz = Dgtd_zxx;
        double Ctd_xyx = Dgtd_yxx;
        double Ctd_xyy = 2*Dgtd_yxy-Dgtd_xyy;
        double Ctd_xyz = t3;
        double Ctd_xzx = Dgtd_zxx;
        double Ctd_xzy = t3;
        double Ctd_xzz = 2*Dgtd_zxz-Dgtd_xzz;
        double Ctd_yxx = 2*Dgtd_xxy-Dgtd_yxx;
        double Ctd_yxy = Dgtd_xyy;
        double Ctd_yxz = t8;
        double Ctd_yyx = Dgtd_xyy;
        double Ctd_yyy = Dgtd_yyy;
        double Ctd_yyz = Dgtd_zyy;
        double Ctd_yzx = t8;
        double Ctd_yzy = Dgtd_zyy;
        double Ctd_yzz = 2*Dgtd_zyz-Dgtd_yzz;
        double Ctd_zxx = 2*Dgtd_xxz-Dgtd_zxx;
        double Ctd_zxy = t13;
        double Ctd_zxz = Dgtd_xzz;
        double Ctd_zyx = t13;
        double Ctd_zyy = 2*Dgtd_yyz-Dgtd_zyy;
        double Ctd_zyz = Dgtd_yzz;
        double Ctd_zzx = Dgtd_xzz;
        double Ctd_zzy = Dgtd_yzz;
        double Ctd_zzz = Dgtd_zzz;
        t8 = gtu_xx*Ctd_xxy+gtu_xy*Ctd_yxy+gtu_xz*Ctd_zxy;
        double t12 = gtu_xx*Ctd_xxz+gtu_xy*Ctd_yxz+gtu_xz*Ctd_zxz;
        double t20 = gtu_xx*Ctd_xyz+gtu_xy*Ctd_yyz+gtu_xz*Ctd_zyz;
        double t32 = gtu_xy*Ctd_xxy+gtu_yy*Ctd_yxy+gtu_yz*Ctd_zxy;
        double t36 = gtu_xy*Ctd_xxz+gtu_yy*Ctd_yxz+gtu_yz*Ctd_zxz;
        double t44 = gtu_xy*Ctd_xyz+gtu_yy*Ctd_yyz+gtu_yz*Ctd_zyz;
        double t56 = gtu_xz*Ctd_xxy+gtu_yz*Ctd_yxy+gtu_zz*Ctd_zxy;
        double t60 = gtu_xz*Ctd_xxz+gtu_yz*Ctd_yxz+gtu_zz*Ctd_zxz;
        double t68 = gtu_xz*Ctd_xyz+gtu_yz*Ctd_yyz+gtu_zz*Ctd_zyz;
        double Ct_xxx = gtu_xx*Ctd_xxx+gtu_xy*Ctd_yxx+gtu_xz*Ctd_zxx;
        double Ct_xxy = t8;
        double Ct_xxz = t12;
        double Ct_xyx = t8;
        double Ct_xyy = gtu_xx*Ctd_xyy+gtu_xy*Ctd_yyy+gtu_xz*Ctd_zyy;
        double Ct_xyz = t20;
        double Ct_xzx = t12;
        double Ct_xzy = t20;
        double Ct_xzz = gtu_xx*Ctd_xzz+gtu_xy*Ctd_yzz+gtu_xz*Ctd_zzz;
        double Ct_yxx = gtu_xy*Ctd_xxx+gtu_yy*Ctd_yxx+gtu_yz*Ctd_zxx;
        double Ct_yxy = t32;
        double Ct_yxz = t36;
        double Ct_yyx = t32;
        double Ct_yyy = gtu_xy*Ctd_xyy+gtu_yy*Ctd_yyy+gtu_yz*Ctd_zyy;
        double Ct_yyz = t44;
        double Ct_yzx = t36;
        double Ct_yzy = t44;
        double Ct_yzz = gtu_xy*Ctd_xzz+gtu_yy*Ctd_yzz+gtu_yz*Ctd_zzz;
        double Ct_zxx = gtu_xz*Ctd_xxx+gtu_yz*Ctd_yxx+gtu_zz*Ctd_zxx;
        double Ct_zxy = t56;
        double Ct_zxz = t60;
        double Ct_zyx = t56;
        double Ct_zyy = gtu_xz*Ctd_xyy+gtu_yz*Ctd_yyy+gtu_zz*Ctd_zyy;
        double Ct_zyz = t68;
        double Ct_zzx = t60;
        double Ct_zzy = t68;
        double Ct_zzz = gtu_xz*Ctd_xzz+gtu_yz*Ctd_yzz+gtu_zz*Ctd_zzz;
        Gamh_x[i3D] = 0.5*(gtu_xx*Ct_xxx+gtu_yy*Ct_xyy+gtu_zz*Ct_xzz)+gtu_xy*Ct_xxy+gtu_xz*Ct_xxz+gtu_yz*Ct_xyz;
        Gamh_y[i3D] = 0.5*(gtu_xx*Ct_yxx+gtu_yy*Ct_yyy+gtu_zz*Ct_yzz)+gtu_xy*Ct_yxy+gtu_xz*Ct_yxz+gtu_yz*Ct_yyz;
        Gamh_z[i3D] = 0.5*(gtu_xx*Ct_zxx+gtu_yy*Ct_zyy+gtu_zz*Ct_zzz)+gtu_xy*Ct_zxy+gtu_xz*Ct_zxz+gtu_yz*Ct_zyz;
      }
    }
  }

  delete[] (dx_gtd11);
  delete[] (dx_gtd12);
  delete[] (dx_gtd13);
  delete[] (dx_gtd22);
  delete[] (dx_gtd23);
  delete[] (dx_gtd33);

  delete[] (dy_gtd11);
  delete[] (dy_gtd12);
  delete[] (dy_gtd13);
  delete[] (dy_gtd22);
  delete[] (dy_gtd23);
  delete[] (dy_gtd33);

  delete[] (dz_gtd11);
  delete[] (dz_gtd12);
  delete[] (dz_gtd13);
  delete[] (dz_gtd22);
  delete[] (dz_gtd23);
  delete[] (dz_gtd33);
}
