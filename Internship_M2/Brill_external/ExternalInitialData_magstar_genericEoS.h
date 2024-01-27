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

static double temp0, temp1;
static double energy_shift;

// min and max values

static double eos_rhomax, eos_rhomin;
static double eos_tempmin, eos_tempmax;
static double eos_yemin, eos_yemax;
static double eos_epsmin, eos_epsmax;
static double eos_pressmin, eos_pressmax;

static double c2p_tempmin;
static double c2p_tempmax;

// table key
// 0 logpress 
// 1 logenergy
// 2 entropy
// 3 munu
// 4 cs2
// 5 dedt
// 6 dpdrhoe
// 7 dpderho
// 8 muhat
// 9 mu_e
// 10 mu_p
// 11 mu_n
// 12 Xa
// 13 Xh
// 14 Xn
// 15 Xp
// 16 Abar
// 17 Zbar
// 18 Gamma
enum eos_var {i_logpress=0, i_logenergy, i_entropy, i_munu, i_cs2, i_dedt,
                i_dpdrhoe, i_dpderho, i_muhat, i_mu_e, i_mu_p, i_mu_n, i_Xa,
                i_Xh, i_Xn, i_Xp, i_Abar, i_Zbar, i_Gamma};

// table data

static int nrho;
static int ntemp;
static int nye;
//
static double* alltables;
static double* epstable;
static double* presstable;
static double* logrho;
static double* logtemp;
static double dlintemp,dlintempi;
static double drholintempi;
static double dlintempyei;
static double drholintempyei;
static double* yes;
static double dtemp, dtempi;
static double drho, drhoi;
static double dye, dyei;
static double drhotempi;
static double drhoyei;
static double dtempyei;
static double drhotempyei;

  static void adm2bssn(double *DYf, double *Ye, double *eps, double *tau, double *rho, double *D,   double *vx,     
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

static void nuc_eos_lookup(double& P, double& eps, double& temp, double& csq, double rho, double ye);
static double ye_beta_equil(double rho0, double temp);
static int be_zbrent(double& x, double x1, double x2, double tol, double rho, double temp);
static double find_min_abs_munu(double& tol, double rho, double temp);
static double be_eqs(double& f, double x, double rho, double temp);
static void check_eostable_input(double xrho,double xtemp,double& xye);
static void nuc_eos_mu(double& xrho,double& xtemp,double xye,double& xmu_e,double& xmu_n,double& xmu_p);
static int nuc_eos_short(double& xrho, double& xtemp, double xye, double& xenr, double& xprs, double& xcs2, double xrfeps);
static void findmus(double* lr,double* lt, double* y, double* ff);
static void findall_short(double* lr,double* lt, double* y, double* ff);
static void findthis(double* lr,double* lt,double* y,double* value, int iv, double& d1, double& d2, double& d3);
static void intp3d(double* x, double* y, double* z, double* f, int kt, double* ft, int* ivs, int nivs, int nx, int ny, int nz, double* xt, double* yt, double* zt, double& d1, double& d2, double& d3);
static void nuc_eos_C_ReadTable(char* nuceos_table_name);
static void nuc_eos_C_UnloadTable();
};



