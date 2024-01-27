
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <sys/stat.h>
#include <cstdlib>

#include "ExternalInitialData.h"

using namespace std;
using namespace Lorene;

#define CUto10KM 1.476887220238973e-01
#define CUtoKM CUto10KM*10.0;
#define CUtoKMBN 1.476887220238973
#define DensityConvFactor 1.619679366684996e-21


#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define CD4(paru, ip2, ip1, im1, im2, dx) ((((-paru[ip2]) + 8.0 * paru[ip1]) + ((-8.0 * paru[im1]) + paru[im2])) / (12.0 * dx))

#define CD6(paru, ip3, ip2, ip1, im1, im2, im3, dx) ((((paru[ip3]) - 9 * paru[ip2] + 45.0 * paru[ip1]) + ((-45.0 * paru[im1]) + 9 * paru[im2] - paru[im3])) / (60.0 * dx))

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
void SAMRAI_External_Data::external_loadData(int vars, ...) {

  va_list args;
  va_start(args, vars);

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
  double* vacuum = va_arg(args, double*);
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

  int rc,i;

  double *x, *y, *z;
  double *nbar, *ener_spec;

  double gamma1, gamma2, kappa1, kappa2, mass1_b, mass2_b;
  double srad1[4], srad2[4], rad1, rad2;

  const double invCUtoKM = 1.0/CUtoKMBN;

  char pfile[256];
  const char *dfile  = "resu.d";

  struct stat buffer;   
  if (stat ("resu.d", &buffer) != 0) {
      std::cout <<"Lorene exception: File resu.d does not exist."<< std::endl;
      std::exit(-1);
  }

  // Allocate storage
  int nd = nx * ny * nz;

  double *x3d = (double *) malloc(nd*sizeof(double));
  double *y3d = (double *) malloc(nd*sizeof(double));
  double *z3d = (double *) malloc(nd*sizeof(double));
  // rescale the coordinates
  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        int i3D = i + nx*(j  + ny*k);
        x3d[i3D] = (init_domain_x + ((i + init_patch_x) - d_ghost_width) * dx[0]) * CUtoKMBN;
        y3d[i3D] = (init_domain_y + ((j + init_patch_y) - d_ghost_width) * dx[1]) * CUtoKMBN;
        z3d[i3D] = (init_domain_z + ((k + init_patch_z) - d_ghost_width) * dx[2]) * CUtoKMBN;
      }
    }
  }     

  // Read Lorene Data
  Bin_NS binary(nd, x3d, y3d, z3d, dfile);

  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        int i3D = i + nx*(j  + ny*k);
        gxx[i3D] = binary.g_xx[i3D];
        gxy[i3D] = binary.g_xy[i3D];
        gxz[i3D] = binary.g_xz[i3D];
        gyy[i3D] = binary.g_yy[i3D];
        gyz[i3D] = binary.g_yz[i3D];
        gzz[i3D] = binary.g_zz[i3D];

        Kxx[i3D] = binary.k_xx[i3D];
        Kxy[i3D] = binary.k_xy[i3D];
        Kxz[i3D] = binary.k_xz[i3D];
        Kyy[i3D] = binary.k_yy[i3D];
        Kyz[i3D] = binary.k_yz[i3D];
        Kzz[i3D] = binary.k_zz[i3D];

        Alpha[i3D]  = binary.nnn[i3D];
        Betaux[i3D] = - binary.beta_x[i3D];
        Betauy[i3D] = - binary.beta_y[i3D];
        Betauz[i3D] = - binary.beta_z[i3D];

        vx[i3D]     = binary.u_euler_x[i3D];  // v^x  (velocity of Eulerian observers)
        vy[i3D]     = binary.u_euler_y[i3D];  // v^y  (velocity of Eulerian observers)
        vz[i3D]     = binary.u_euler_z[i3D];  // v^z  (velocity of Eulerian observers)
      }
    }
  }

  nbar      = binary.nbar;       // baryon density, kg/m^3
  ener_spec = binary.ener_spec;  // specific internal energy, c^2
  
  gamma1 = binary.gamma_poly1;
  gamma2 = binary.gamma_poly2;
  kappa1 = binary.kappa_poly1;
  kappa2 = binary.kappa_poly2;
  mass1_b  = binary.mass1_b;
  mass2_b  = binary.mass2_b;

  double dist    = binary.dist        * invCUtoKM;

  // Read in the four radii for the stars
  srad1[0]  = binary.rad1_x_comp;
  srad1[1]  = binary.rad1_y;
  srad1[2]  = binary.rad1_z;
  srad1[3]  = binary.rad1_x_opp;

  srad2[0]  = binary.rad2_x_comp;
  srad2[1]  = binary.rad2_y;
  srad2[2]  = binary.rad2_z;
  srad2[3]  = binary.rad2_x_opp;

  // Find the maximum radius.  This is used to add B field 
  rad1 = 0.0;
  rad2 = 0.0;
  for (int ii = 0; ii < 4; ii++) 
    if (srad1[ii] > rad1) rad1 = srad1[ii];
  for (int ii = 0; ii < 4; ii++) 
    if (srad2[ii] > rad2) rad2 = srad2[ii];

  double r1 = rad1 * invCUtoKM;
  double r2 = rad2 * invCUtoKM;


  if (fabs(gamma1-gamma2) < 1.0e-8) {
    *Gamma = gamma1;
  }
  else {
    cerr << "gamma1 != gamma2" << endl;
    cerr << "   gamma1" << gamma1 << endl;
    cerr << "   gamma2" << gamma2 << endl;
    return;
  }

  // Some unit conversions and variable transformations 
  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        int i3D = i + nx*(j  + ny*k);
        Kxx[i3D] *= CUtoKMBN;
        Kxy[i3D] *= CUtoKMBN;
        Kxz[i3D] *= CUtoKMBN;
        Kyy[i3D] *= CUtoKMBN;
        Kyz[i3D] *= CUtoKMBN;
        Kzz[i3D] *= CUtoKMBN;

        rho[i3D] = nbar[i3D]*DensityConvFactor;
        eps[i3D] = ener_spec[i3D];

        if (rho[i3D] < (*vacuum)) rho[i3D] = (*vacuum);
        if (eps[i3D]   < (*vacuum)) eps[i3D] = (*vacuum);

        double p_cold, eps_cold, selected_gamma;
        if (rho[i3D] > (*rho_2)) {
          selected_gamma = (*gamma_3);
          p_cold = (*K_3) * (pow(rho[i3D], (*gamma_3)));
          eps_cold = (*a_3) + (*K_3) / ((*gamma_3) - 1.0) * (pow(rho[i3D], ((*gamma_3) - 1.0)));
        } else {
          if (!(*rho_1 > rho[i3D])) {
            selected_gamma = (*gamma_2);
            p_cold = (*K_2) * (pow(rho[i3D], (*gamma_2)));
            eps_cold = (*a_2) + (*K_2) / ((*gamma_2) - 1.0) * (pow(rho[i3D], ((*gamma_2) - 1.0)));
          } else {
            if (!(*rho_0 > rho[i3D])) {
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

        P[i3D] = p_cold + ((*Gamma) - 1) * rho[i3D] * (eps[i3D] - eps_cold);

        // Apply floor 
        if (P[i3D]   < (*vacuum)*(*vacuum)) P[i3D] = (*vacuum)*(*vacuum);

        double h = rho[i3D] * (1.0 + eps[i3D]) + P[i3D];
        sqcs[i3D] = fabs(((*Gamma) * P[i3D] + (selected_gamma - (*Gamma)) * p_cold) / h);

        // Explicitly set the B-field to zero.
        Bx[i3D] = 0.0;
        By[i3D] = 0.0;
        Bz[i3D] = 0.0;

        // Explicitly set the Shift-field to zero.    
        Betaux[i3D] = 0.0;
        Betauy[i3D] = 0.0;
        Betauz[i3D] = 0.0;
        
        theta[i3D] = 0;

      }
    }
  }

  free(x3d);
  free(y3d);
  free(z3d);

  magneticField(P, Sx, Sy, Sz, Bx, By, Bz, gxx, gxy, gxz, gyy, gyz, gzz, *tov_angle, *tov_Asize, *tov_facp, *vacuum, r1, r2, dist, init_domain_x, init_domain_y, init_domain_z, init_patch_x, init_patch_y, init_patch_z, d_ghost_width, dx, nx, ny, nz);

  adm2bssn(eps, tau, rho, D, vx, 
           vy, vz, Sx, Sy, Sz, 
           P, Bx, By, Bz, gxx, 
           gxy, gxz, gyy, gyz, gzz,
           Kxx, Kxy, Kxz, Kyy, Kyz, 
           Kzz, Alpha, chi, trK, theta, 
           Gamh_x, Gamh_y, Gamh_z, *Gamma, Betaux, 
           Betauy, Betauz, dx, nx, ny, 
           nz);
}

void SAMRAI_External_Data::magneticField(double *P,   double *Sx,     double *Sy,     double *Sz, double *Bx,     double *By,     double *Bz, 
              double *gxx,   double *gxy,    double *gxz,
              double *gyy,   double *gyz,    double *gzz,
              double tov_angle, double tov_Asize, double tov_facp, double vacuum, double rad1, double rad2, double dist, double init_domain_x, double init_domain_y, double init_domain_z, int init_patch_x, int init_patch_y, int init_patch_z, int d_ghost_width, const double* dxA, int nx,        int ny,         int nz) {

  double dx = dxA[0];
  double dy = dxA[1];
  double dz = dxA[2];

  /* location of star 1 */
  double xs1 = 0.5*dist;
  double ys1 = 0.0;
  double zs1 = 0.0;
  /* location of star 2 */
  double xs2 = -0.5*dist;
  double ys2 = 0.0;
  double zs2 = 0.0;

  double t1 = pow(1.4*rad1,2);
  double t2 = pow(1.4*rad2,2);

  double r1sq, r2sq, x_fd, y_fd, z_fd;
  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        int i3D = i + nx*(j  + ny*k);

        // convert position from cactus units to Lorene units (10Km)
        x_fd = (init_domain_x + ((i + init_patch_x) - d_ghost_width) * dxA[0]);
        y_fd = (init_domain_y + ((j + init_patch_y) - d_ghost_width) * dxA[1]);
        z_fd = (init_domain_z + ((k + init_patch_z) - d_ghost_width) * dxA[2]);


        r1sq = (x_fd-xs1)*(x_fd-xs1) + (y_fd-ys1)*(y_fd-ys1) + (z_fd-zs1)*(z_fd-zs1);
        r2sq = (x_fd-xs2)*(x_fd-xs2) + (y_fd-ys2)*(y_fd-ys2) + (z_fd-zs2)*(z_fd-zs2);

        if (r1sq <= t1) {
          Sx[i3D]  = -(y_fd-ys1)*tov_Asize*max(P[i3D]-tov_facp*vacuum, 0.0);
          Sy[i3D]  =  (x_fd-xs1)*tov_Asize*max(P[i3D]-tov_facp*vacuum, 0.0);
          Sz[i3D]  = 0.0;
        }
        if (r2sq <= t2) {
            Sx[i3D]  = -(y_fd-ys2)*tov_Asize*max(P[i3D]-tov_facp*vacuum, 0.0);
            Sy[i3D]  =  (x_fd-xs2)*tov_Asize*max(P[i3D]-tov_facp*vacuum, 0.0);
            Sz[i3D]  = 0.0;
        }
        if (r1sq > t1 && r2sq > t2) {
          Sx[i3D]  = 0.0;
          Sy[i3D]  = 0.0;
          Sz[i3D]  = 0.0;
        }


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

void SAMRAI_External_Data::adm2bssn(double *eps, double *tau, double *rho, double *D,   double *vx,     
                                    double *vy,     double *vz,   double *Sx,     double *Sy,  double *Sz,
                                    double *P,     double *Bx,     double *By,     double *Bz, double *gxx,   
                                    double *gxy,    double *gxz, double *gyy,   double *gyz,    double *gzz,
                                    double *Kxx,   double *Kxy,    double *Kxz, double *Kyy,   double *Kyz,    
                                    double *Kzz, double *Alpha, double *chi, double *trK, double *theta,
                                    double *Gamh_x,     double *Gamh_y,     double *Gamh_z, double Gamma, double *Betaux,
                                    double *Betauy, double *Betauz, const double* dxA, int nx,        int ny,         
                                    int nz) {

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

  double *dx_gtd11 = (double *) malloc(nd*sizeof(double));
  double *dx_gtd12 = (double *) malloc(nd*sizeof(double));
  double *dx_gtd13 = (double *) malloc(nd*sizeof(double));
  double *dx_gtd22 = (double *) malloc(nd*sizeof(double));
  double *dx_gtd23 = (double *) malloc(nd*sizeof(double));
  double *dx_gtd33 = (double *) malloc(nd*sizeof(double));

  double *dy_gtd11 = (double *) malloc(nd*sizeof(double));
  double *dy_gtd12 = (double *) malloc(nd*sizeof(double));
  double *dy_gtd13 = (double *) malloc(nd*sizeof(double));
  double *dy_gtd22 = (double *) malloc(nd*sizeof(double));
  double *dy_gtd23 = (double *) malloc(nd*sizeof(double));
  double *dy_gtd33 = (double *) malloc(nd*sizeof(double));

  double *dz_gtd11 = (double *) malloc(nd*sizeof(double));
  double *dz_gtd12 = (double *) malloc(nd*sizeof(double));
  double *dz_gtd13 = (double *) malloc(nd*sizeof(double));
  double *dz_gtd22 = (double *) malloc(nd*sizeof(double));
  double *dz_gtd23 = (double *) malloc(nd*sizeof(double));
  double *dz_gtd33 = (double *) malloc(nd*sizeof(double));

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

  free(dx_gtd11);
  free(dx_gtd12);
  free(dx_gtd13);
  free(dx_gtd22);
  free(dx_gtd23);
  free(dx_gtd33);

  free(dy_gtd11);
  free(dy_gtd12);
  free(dy_gtd13);
  free(dy_gtd22);
  free(dy_gtd23);
  free(dy_gtd33);

  free(dz_gtd11);
  free(dz_gtd12);
  free(dz_gtd13);
  free(dz_gtd22);
  free(dz_gtd23);
  free(dz_gtd33);
}
