/*-------------------------------------------------------------------------
 *
 *   DWN (08/06/13): Some updates for nuclear EOS
 *   SLL (06/07/13): changing for realistic EOS
 *                   by looking at lorene_star_rot.cc for what to change
 *-------------------------------------------------------------------------*/
#include "et_rot_mag.h"
#include "eos.h"
#include "nbr_spx.h"
#include "unites.h"
#include "metric.h"


using namespace std;
using namespace Lorene;

class SAMRAI_External_Data
{

public:

   static void external_loadData(int vars, ...);

private:

  static void adm2bssn(double *eps, double *tau, double *rho, double *D,   double *vx,     
                        double *vy,     double *vz,   double *Sx,     double *Sy,     double *Sz,
                double *P,     double *Bx,     double *By,     double *Bz, double *gxx,   
                double *gxy,    double *gxz, double *gyy,   double *gyz,    double *gzz,
                double *Kxx,   double *Kxy,    double *Kxz, double *Kyy,   double *Kyz,    
                double *Kzz, double *Alpha, double *chi, double *trK, double *theta,
                double *Gamh_x,     double *Gamh_y,     double *Gamh_z, double Gamma, double *Betaux,
                double *Betauy, double *Betauz, const double* dxA, int nx,        int ny,         int nz);

  static void magneticField(double *P,   double *Sx,     double *Sy,     double *Sz, double *Bx,     double *By,     double *Bz, 
              double *gxx,   double *gxy,    double *gxz,
              double *gyy,   double *gyz,    double *gzz,
              double tov_angle, double tov_Asize, double tov_facp, double vacuum, double init_domain_x, double init_domain_y, double init_domain_z, int init_patch_x, int init_patch_y, int init_patch_z, int d_ghost_width, const double* dxA, int nx,        int ny,         int nz);

};



