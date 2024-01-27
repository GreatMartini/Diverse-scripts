
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <sys/stat.h>
#include <cstdlib>

#include "ExternalInitialData.h"
#include "unites.h"


using namespace std;
using namespace Lorene;
using namespace Unites;

const int EXTRAP_POLY_NUM_POINTS = 8;
const double POLY_EXTRAP_DR = 0.10;

typedef struct {
   int i, j, k;
   double x, y, z;
} AH_Node_Point;



#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define CD4(paru, ip2, ip1, im1, im2, dx) ((((-paru[ip2]) + 8.0 * paru[ip1]) + ((-8.0 * paru[im1]) + paru[im2])) / (12.0 * dx))

#define CD6(paru, ip3, ip2, ip1, im1, im2, im3, dx) ((((paru[ip3]) - 9 * paru[ip2] + 45.0 * paru[ip1]) + ((-45.0 * paru[im1]) + 9 * paru[im2] - paru[im3])) / (60.0 * dx))

/*double *eps,    double *tau,    double *rho,    double *D,   
  double *vx,     double *vy,     double *vz,   
  double *Sx,     double *Sy,     double *Sz,
  double *P,     
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

  double* eps3d = va_arg(args, double *);
  double* tau = va_arg(args, double *);
  double* rho3d = va_arg(args, double *);
  double* D = va_arg(args, double *);
  double* vx3d = va_arg(args, double *);
  double* vy3d = va_arg(args, double *);
  double* vz3d = va_arg(args, double *);
  double* Sx = va_arg(args, double *);
  double* Sy = va_arg(args, double *);
  double* Sz = va_arg(args, double *);
  double* P = va_arg(args, double *);
  double* sqcs = va_arg(args, double *);
  double* Bx = va_arg(args, double *);
  double* By = va_arg(args, double *);
  double* Bz = va_arg(args, double *);
  double* gxx3d = va_arg(args, double *);
  double* gxy3d = va_arg(args, double *);
  double* gxz3d = va_arg(args, double *);
  double* gyy3d = va_arg(args, double *);
  double* gyz3d = va_arg(args, double *);
  double* gzz3d = va_arg(args, double *);
  double* Kxx3d = va_arg(args, double *);
  double* Kxy3d = va_arg(args, double *);
  double* Kxz3d = va_arg(args, double *);
  double* Kyy3d = va_arg(args, double *);
  double* Kyz3d = va_arg(args, double *);
  double* Kzz3d = va_arg(args, double *);
  double* Alpha3d = va_arg(args, double *);
  double* chi = va_arg(args, double *);
  double* trK = va_arg(args, double *);
  double* theta = va_arg(args, double *);
  double* Gamh_x = va_arg(args, double *);
  double* Gamh_y = va_arg(args, double *);
  double* Gamh_z = va_arg(args, double *);
  double* Betaux3d = va_arg(args, double *);
  double* Betauy3d = va_arg(args, double *);
  double* Betauz3d = va_arg(args, double *);
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

  double sigma = 1.0;
  double amp = 0.1;
  AH_Node_Point ahpts[2000];

// Read in the data file.
  const char* dataFilename = "bin.dat";
  struct stat buffer;   
  if (stat (dataFilename, &buffer) != 0) {
      std::cout <<"Lorene exception: File resu.d does not exist."<< std::endl;
      std::exit(-1);
  }

  FILE* filePointer = fopen(dataFilename, "r") ;
  Mg3d mg_ns(filePointer) ;
  Map_et mp_ns(mg_ns, filePointer) ;
  Eos* peos = Eos::eos_from_file(filePointer) ;
  Mg3d mg_bh(filePointer) ;
  Map_af mp_bh(mg_bh, filePointer) ;
  Bin_ns_bh theBinary(mp_ns, *peos, mp_bh, filePointer) ;    
  fclose(filePointer) ;
  
  // The Following lines are required to initialise various objects
  // insde the binary bh ns object.
  theBinary.set_ns().update_metric( theBinary.get_bh() ) ;
  theBinary.set_bh().fait_n_comp( theBinary.get_ns() ) ;
  theBinary.set_bh().fait_psi_comp( theBinary.get_ns() ) ;
  theBinary.set_bh().fait_taij_auto( ) ;
  theBinary.set_ns().update_metric_der_comp( theBinary.get_bh() ) ;
  theBinary.set_bh().update_metric (theBinary.get_ns()) ;
  theBinary.fait_tkij() ;
  theBinary.set_ns().equation_of_state() ;
  theBinary.set_ns().kinematics(theBinary.get_omega(), theBinary.get_x_axe()) ;
  theBinary.set_ns().fait_d_psi() ;
  theBinary.set_ns().hydro_euler() ;

  // Stuff to store important quantities in, temporarily.
  // Mainly only necessary for the storage, other than the x,y,z arryas.
  Et_bin_nsbh ns = theBinary.get_ns();
  Bhole bh = theBinary.get_bh();

  // Get the tensors corresponding to the fluid variables.
  Eos_poly* castToPolytrope = dynamic_cast<Eos_poly*>(peos);
  Cmp tensorRho = (*castToPolytrope).get_m_0() * ns.get_nbar()();
  Cmp tensorEpsilon = ns.get_ener()() - tensorRho;
  // ^ Not specific internal energy (epsilon), but rho*eps.
  Cmp tensorVx = ns.get_u_euler()(0);
  Cmp tensorVy = ns.get_u_euler()(1);
  Cmp tensorVz = ns.get_u_euler()(2);

  // Spacetime variables, for the ns and the bh.
  Cmp lapseNS = ns.get_n_auto()();
  Cmp lapseBH = bh.get_n_auto()();
  Cmp shiftxNS = ns.get_shift_auto()(0);
  Cmp shiftyNS = ns.get_shift_auto()(1);
  Cmp shiftzNS = ns.get_shift_auto()(2);
  Cmp shiftxBH = bh.get_shift_auto()(0);
  Cmp shiftyBH = bh.get_shift_auto()(1);
  Cmp shiftzBH = bh.get_shift_auto()(2);
  // Psi is the conformal factor.
  Cmp psiNS = ns.get_confpsi_auto()();
  // get_psi_auto would return the velocity potential.
  Cmp psiBH = bh.get_psi_auto()();

  // The extrinsic curvature, a bunch of necessary gibberish.
  // Honestly, I don't know why this is not handled automatically.
  //
  // I am not at all sure that the following lines are required to
  //  properly extract the extrinsic curvature, but where the extrinsic
  //  curvature is stored in contravariant form, I am just unwilling to
  //  take the time to try to decipher the following code. For that reason
  //  what follows is simply a copy of what I have seen in other codes.
  Tenseur_sym aNS(ns.get_taij_auto());
  aNS.set_std_base();
  aNS.dec2_dzpuis();
  Tenseur_sym aBH(bh.get_taij_auto());
  aBH.set_std_base();
  aBH.dec2_dzpuis();
  // This garbage makes sense, in one frame of thinking. We need to
  //  pull out the old tensor objects, because the new tensor objects
  //  don't have get_value implemented.
  Valeur axxNS = (aNS(0,0)).va;
  Valeur axxBH = (aBH(0,0)).va;
  Valeur axyNS = (aNS(0,1)).va;
  Valeur axyBH = (aBH(0,1)).va;
  Valeur axzNS = (aNS(0,2)).va;
  Valeur axzBH = (aBH(0,2)).va;
  Valeur ayyNS = (aNS(1,1)).va;
  Valeur ayyBH = (aBH(1,1)).va;
  Valeur ayzNS = (aNS(1,2)).va;
  Valeur ayzBH = (aBH(1,2)).va;
  Valeur azzNS = (aNS(2,2)).va;
  Valeur azzBH = (aBH(2,2)).va;
  // Now for some reason, we also need to set the spectral basis for
  //  these tensors. Who knows.
  axxNS.set_base( (ns.get_taij_auto()(0,0)).va.base );
  axxBH.set_base( (bh.get_taij_auto()(0,0)).va.base );
  axyNS.set_base( (ns.get_taij_auto()(0,1)).va.base );
  axyBH.set_base( (bh.get_taij_auto()(0,1)).va.base );
  axzNS.set_base( (ns.get_taij_auto()(0,2)).va.base );
  axzBH.set_base( (bh.get_taij_auto()(0,2)).va.base );
  ayyNS.set_base( (ns.get_taij_auto()(1,1)).va.base );
  ayyBH.set_base( (bh.get_taij_auto()(1,1)).va.base );
  ayzNS.set_base( (ns.get_taij_auto()(1,2)).va.base );
  ayzBH.set_base( (bh.get_taij_auto()(1,2)).va.base );
  azzNS.set_base( (ns.get_taij_auto()(2,2)).va.base );
  azzBH.set_base( (bh.get_taij_auto()(2,2)).va.base );
  // Recall that the above extrinsic curvature is in contravariant form.
  // Now for some more black magic: compute the conefficients!
  axxNS.coef();
  axxBH.coef();
  axyNS.coef();
  axyBH.coef();
  axzNS.coef();
  axzBH.coef();
  ayyNS.coef();
  ayyBH.coef();
  ayzNS.coef();
  ayzBH.coef();
  azzNS.coef();
  azzBH.coef();

  double CUtoLUlength = g_si*msol_si/pow(c_si,2)/r_unit;
  double CUto10km = CUtoLUlength;
  double LUtoCUrho = pow(c_si,-6)*pow(msol_si,2)*pow(g_si,3)*rho_unit;

  double rho, eps, velx, vely, velz, alp, betax, betay, betaz;
  double psi, gxx, gxy, gxz, gyy, gyz, gzz;
  double      kxx, kxy, kxz, kyy, kyz, kzz;
  double      axx, axy, axz, ayy, ayz, azz;

  double bh_radius = bh.get_rayon();
  double bh_origin_x = bh.get_mp().get_ori_x();
  double bh_origin_y = bh.get_mp().get_ori_y();
  double bh_origin_z = bh.get_mp().get_ori_z();

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
        x3d[i3D] = (init_domain_x + ((i + init_patch_x) - d_ghost_width) * dx[0]);
        y3d[i3D] = (init_domain_y + ((j + init_patch_y) - d_ghost_width) * dx[1]);
        z3d[i3D] = (init_domain_z + ((k + init_patch_z) - d_ghost_width) * dx[2]);
      }
    }
  }  

  int ahcount = 0;
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        int i3D = i + nx*(j  + ny*k);
        double xLoop = x3d[i3D]*CUto10km;
        double yLoop = y3d[i3D]*CUto10km;
        double zLoop = z3d[i3D]*CUto10km;
        double radBH = pow(xLoop - bh_origin_x,2.0) 
                     + pow(yLoop - bh_origin_y,2.0)
                     + pow(zLoop - bh_origin_z,2.0);
        radBH = sqrt(radBH);
        double chi;
        if (radBH < bh_radius) {
          radBH = radBH/bh_radius*acos(-1.0e0)/2.0e0;
          radBH = -tan(radBH)*tan(radBH)/sigma/sigma;
          chi = 1.0e0 - amp*exp(radBH);
          if (ahcount == 2000) {
            cerr << "Ran out of memory in ahpts" << endl;
            exit(2);
          }
          ahpts[ahcount].i = i;
          ahpts[ahcount].j = j;
          ahpts[ahcount].k = k;
          ahpts[ahcount].x = xLoop;
          ahpts[ahcount].y = yLoop;
          ahpts[ahcount].z = zLoop;
          ahcount++;
        } else {
          chi = 1.0e0;
        }

        // Local coordinates for the BH and NS spectral domain.
        double rNS, thetaNS, phiNS;
        double rBH, thetaBH, phiBH;
        mp_ns.convert_absolute(xLoop,yLoop,zLoop, rNS, thetaNS, phiNS);
        mp_bh.convert_absolute(xLoop,yLoop,zLoop, rBH, thetaBH, phiBH);

        // Fluid variables alone.
        rho =     tensorRho.val_point(rNS,thetaNS,phiNS);
        if(rho>0) {
          eps = tensorEpsilon.val_point(rNS,thetaNS,phiNS)/rho;
        } else {
          eps = 0.0;
        }

        // Now converted to normal specific internal energy.
        velx = tensorVx.val_point(rNS,thetaNS,phiNS);
        vely = tensorVy.val_point(rNS,thetaNS,phiNS);
        velz = tensorVz.val_point(rNS,thetaNS,phiNS);

        // Lapse and shift.
        alp = lapseNS.val_point(rNS,thetaNS,phiNS) +
                 lapseBH.val_point(rBH,thetaBH,phiBH);
        betax = shiftxNS.val_point(rNS,thetaNS,phiNS) +
             shiftxBH.val_point(rBH,thetaBH,phiBH);
        betay = shiftyNS.val_point(rNS,thetaNS,phiNS) +
             shiftyBH.val_point(rBH,thetaBH,phiBH);
        betaz = shiftzNS.val_point(rNS,thetaNS,phiNS) +
             shiftzBH.val_point(rBH,thetaBH,phiBH);

        // Conformal factor and metric quantities.
        psi = psiNS.val_point(rNS,thetaNS,phiNS) +
                 psiBH.val_point(rBH,thetaBH,phiBH);
        gxx = pow(psi,4);
        gxy = 0.0e0;
        gxz = 0.0e0;
        gyy = pow(psi,4);
        gyz = 0.0e0;
        gzz = pow(psi,4);
      
        // Extrinsic curvature. We need to retrieve the values in another
        //  bizarre manner. For this we need a couple auxiliary variables.
        int indexNS, indexBH; // Indexes which domain the coordinates lay in.
        double radialNS, radialBH; // The radial coordinate in that domain.
        // Now we initialise them using the object's methods.
        mp_ns.val_lx(rNS,thetaNS,phiNS,indexNS,radialNS);
        mp_bh.val_lx(rBH,thetaBH,phiBH,indexBH,radialBH);
        // Now we use the special value functions to acquire what we need.
        axx = axxNS.c_cf->val_point_asymy(indexNS,radialNS,thetaNS,phiNS) +
             axxBH.c_cf->val_point_asymy(indexBH,radialBH,thetaBH,phiBH);
        axy = axyNS.c_cf->val_point_symy( indexNS,radialNS,thetaNS,phiNS) +
             axyBH.c_cf->val_point_symy( indexBH,radialBH,thetaBH,phiBH);
        axz = axzNS.c_cf->val_point_asymy(indexNS,radialNS,thetaNS,phiNS) +
             axzBH.c_cf->val_point_asymy(indexBH,radialBH,thetaBH,phiBH);
        ayy = ayyNS.c_cf->val_point_asymy(indexNS,radialNS,thetaNS,phiNS) +
             ayyBH.c_cf->val_point_asymy(indexBH,radialBH,thetaBH,phiBH);
        ayz = ayzNS.c_cf->val_point_symy( indexNS,radialNS,thetaNS,phiNS) +
             ayzBH.c_cf->val_point_symy( indexBH,radialBH,thetaBH,phiBH);
        azz = azzNS.c_cf->val_point_asymy(indexNS,radialNS,thetaNS,phiNS) +
             azzBH.c_cf->val_point_asymy(indexBH,radialBH,thetaBH,phiBH);

        // Calculate the extrinsic curvature from A^{ij}.  This is the
        // conformal extrinsic curvature with indices raised or not.
        // Since the data are conformally flat, indices are raised and
        // lowered with the identity.
        kxx = axx / (2.0*alp);
        kxy = axy / (2.0*alp);
        kxz = axz / (2.0*alp);
        kyy = ayy / (2.0*alp);
        kyz = ayz / (2.0*alp);
        kzz = azz / (2.0*alp);

        // To get the ADM (non-conformal) extrinsic curvature, we multiply
        // by \Psi^4.  gxx=\Psi^4, so we multiply by this factor.
        kxx *= gxx;
        kxy *= gxx;
        kxz *= gxx;
        kyy *= gxx;
        kyz *= gxx;
        kzz *= gxx;

        // Next we will convert units using the conversion factors defined
        //   outside the top level loop.
        // The units of the extrinsic curvature are 1/L
        rho *= LUtoCUrho;
        kxx *= CUtoLUlength;
        kxy *= CUtoLUlength;
        kxz *= CUtoLUlength;
        kyy *= CUtoLUlength;
        kyz *= CUtoLUlength;
        kzz *= CUtoLUlength;

        rho3d[i3D] = rho;
        vx3d[i3D]  = velx;
        vy3d[i3D]  = vely;
        vz3d[i3D]  = velz;
        eps3d[i3D] = eps;

        gxx3d[i3D] = gxx;
        gxy3d[i3D] = gxy;
        gxz3d[i3D] = gxz;
        gyy3d[i3D] = gyy;
        gyz3d[i3D] = gyz;
        gzz3d[i3D] = gzz;
        Kxx3d[i3D] = kxx;
        Kxy3d[i3D] = kxy;
        Kxz3d[i3D] = kxz;
        Kyy3d[i3D] = kyy;
        Kyz3d[i3D] = kyz;
        Kzz3d[i3D] = kzz;
        Alpha3d[i3D] = alp;
        Betaux3d[i3D] = betax;
        Betauy3d[i3D] = betay;
        Betauz3d[i3D] = betaz;

      }
    }
  }
   
  /* Check to see if points inside the horizon need to be filled by
   * extrapolation
   */
  if (ahcount > 0) {
    cerr << "Processing " << ahcount << " points inside AH" << endl;
    cerr << "radius of the AH is " << bh_radius << endl;
    double Kxx_extrap[EXTRAP_POLY_NUM_POINTS];
    double Kxy_extrap[EXTRAP_POLY_NUM_POINTS];
    double Kxz_extrap[EXTRAP_POLY_NUM_POINTS];
    double Kyy_extrap[EXTRAP_POLY_NUM_POINTS];
    double Kyz_extrap[EXTRAP_POLY_NUM_POINTS];
    double Kzz_extrap[EXTRAP_POLY_NUM_POINTS];
    double psi_extrap[EXTRAP_POLY_NUM_POINTS];
    double Alpha_extrap[EXTRAP_POLY_NUM_POINTS];
    double Betaux_extrap[EXTRAP_POLY_NUM_POINTS];
    double Betauy_extrap[EXTRAP_POLY_NUM_POINTS];
    double Betauz_extrap[EXTRAP_POLY_NUM_POINTS];
    double r_pts[EXTRAP_POLY_NUM_POINTS];
    double Kxx_new, Kxy_new, Kxz_new, Kyy_new, Kyz_new, Kzz_new;
    double Alpha_new, Betaux_new, Betauy_new, Betauz_new, psi_new;

    for(int ii = 0; ii < ahcount; ii++) {
      /* coords centered on BH */
//        cerr << "... AH point " << ii;
      double cx = ahpts[ii].x - bh_origin_x;
      double cy = ahpts[ii].y - bh_origin_y;
      double cz = ahpts[ii].z - bh_origin_z;
      double cr = sqrt(cx*cx + cy*cy + cz*cz);
//        cerr << "... at (" << cx << ", " << cy << ", " << cz << ")." << endl;
//        cerr << "... at" << cr << endl;
       
      double phi   = atan2(cy,cx);
      double theta = atan2(sqrt(cx*cx + cy*cy),cz);

      for (int kk = 0; kk < EXTRAP_POLY_NUM_POINTS; kk++) {
        double rr = bh_radius + POLY_EXTRAP_DR * ((double)kk + 1.0e-8);
        double xLoop = rr*sin(theta)*cos(phi) + bh_origin_x;
        double yLoop = rr*sin(theta)*sin(phi) + bh_origin_y;
        double zLoop = rr*cos(theta) + bh_origin_z;
        r_pts[kk] = rr;

        // Local coordinates for the BH and NS spectral domain.
        double rNS, thetaNS, phiNS;
        double rBH, thetaBH, phiBH;
        mp_ns.convert_absolute(xLoop,yLoop,zLoop, rNS, thetaNS, phiNS);
              mp_bh.convert_absolute(xLoop,yLoop,zLoop, rBH, thetaBH, phiBH);

        // Lapse and shift.
        alp = lapseNS.val_point(rNS,thetaNS,phiNS) +
                       lapseBH.val_point(rBH,thetaBH,phiBH);
        betax = shiftxNS.val_point(rNS,thetaNS,phiNS) +
                   shiftxBH.val_point(rBH,thetaBH,phiBH);
        betay = shiftyNS.val_point(rNS,thetaNS,phiNS) +
                   shiftyBH.val_point(rBH,thetaBH,phiBH);
        betaz = shiftzNS.val_point(rNS,thetaNS,phiNS) +
                   shiftzBH.val_point(rBH,thetaBH,phiBH);

        // Conformal factor and metric quantities.
        psi = psiNS.val_point(rNS,thetaNS,phiNS) +
                 psiBH.val_point(rBH,thetaBH,phiBH);
          
        // Extrinsic curvature. We need to retrieve the values in another
        //  bizarre manner. For this we need a couple auxiliary variables.
        int indexNS, indexBH; // Indexes which domain the coordinates lay in.
        double radialNS, radialBH; // The radial coordinate in that domain.
        // Now we initialise them using the object's methods.
        mp_ns.val_lx(rNS,thetaNS,phiNS,indexNS,radialNS);
        mp_bh.val_lx(rBH,thetaBH,phiBH,indexBH,radialBH);
        // Now we use the special value functions to acquire what we need.
        axx = axxNS.c_cf->val_point_asymy(indexNS,radialNS,thetaNS,phiNS) +
             axxBH.c_cf->val_point_asymy(indexBH,radialBH,thetaBH,phiBH);
        axy = axyNS.c_cf->val_point_symy( indexNS,radialNS,thetaNS,phiNS) +
             axyBH.c_cf->val_point_symy( indexBH,radialBH,thetaBH,phiBH);
        axz = axzNS.c_cf->val_point_asymy(indexNS,radialNS,thetaNS,phiNS) +
             axzBH.c_cf->val_point_asymy(indexBH,radialBH,thetaBH,phiBH);
        ayy = ayyNS.c_cf->val_point_asymy(indexNS,radialNS,thetaNS,phiNS) +
             ayyBH.c_cf->val_point_asymy(indexBH,radialBH,thetaBH,phiBH);
        ayz = ayzNS.c_cf->val_point_symy( indexNS,radialNS,thetaNS,phiNS) +
             ayzBH.c_cf->val_point_symy( indexBH,radialBH,thetaBH,phiBH);
        azz = azzNS.c_cf->val_point_asymy(indexNS,radialNS,thetaNS,phiNS) +
             azzBH.c_cf->val_point_asymy(indexBH,radialBH,thetaBH,phiBH);

        // Calculate the extrinsic curvature from A^{ij}.  This is the
        // conformal extrinsic curvature with indices raised or not.
        // Since the data are conformally flat, indices are raised and
        // lowered with the identity.
        kxx = axx / (2.0*alp);
        kxy = axy / (2.0*alp);
        kxz = axz / (2.0*alp);
        kyy = ayy / (2.0*alp);
        kyz = ayz / (2.0*alp);
        kzz = azz / (2.0*alp);

        // To get the ADM (non-conformal) extrinsic curvature, we multiply
        // by \Psi^4.  gxx=\Psi^4, so we multiply by this factor.
        gxx = pow(psi,4);
        kxx *= gxx;
        kxy *= gxx;
        kxz *= gxx;
        kyy *= gxx;
        kyz *= gxx;
        kzz *= gxx;

        Kxx_extrap[kk] = kxx;
        Kxy_extrap[kk] = kxy;
        Kxz_extrap[kk] = kxz;
        Kyy_extrap[kk] = kyy;
        Kyz_extrap[kk] = kyz;
        Kzz_extrap[kk] = kzz;
        Alpha_extrap[kk] = alp;
        Betaux_extrap[kk] = betax;
        Betauy_extrap[kk] = betay;
        Betauz_extrap[kk] = betaz;
        psi_extrap[kk] = psi;
      }
      /* Do the extrapolation */  
      double dydummy;
//        for (int jj = 0; jj < EXTRAP_POLY_NUM_POINTS; jj++) {
//          cerr << "Point " << jj << " at r = " << r_pts[jj]
//               << " psi = " << psi_extrap[jj] << " (AH at "
//               << bh_radius << ")" << endl;
//        }
      polint(r_pts, psi_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &psi_new,    &dydummy);
//        cerr << "At r = " << cr << " psi_new = " << psi_new << endl;

      polint(r_pts, Kxx_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Kxx_new,    &dydummy);
      polint(r_pts, Kxy_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Kxy_new,    &dydummy);
      polint(r_pts, Kxz_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Kxz_new,    &dydummy);
      polint(r_pts, Kyy_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Kyy_new,    &dydummy);
      polint(r_pts, Kyz_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Kyz_new,    &dydummy);
      polint(r_pts, Kzz_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Kzz_new,    &dydummy);

      polint(r_pts, Alpha_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Alpha_new,    &dydummy);
      polint(r_pts, Betaux_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Betaux_new,    &dydummy);
      polint(r_pts, Betauy_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Betauy_new,    &dydummy);
      polint(r_pts, Betauz_extrap, EXTRAP_POLY_NUM_POINTS, cr, 
                            &Betauz_new,    &dydummy);

      int i3D = ahpts[ii].i + nx*(ahpts[ii].j  + ny*ahpts[ii].k);
      gxx3d[i3D] = pow(psi_new,4);
      gxy3d[i3D] = 0.0;
      gxz3d[i3D] = 0.0;
      gyy3d[i3D] = pow(psi_new,4);
      gyz3d[i3D] = 0.0;
      gzz3d[i3D] = pow(psi_new,4);
      Kxx3d[i3D] = Kxx_new * CUtoLUlength;
      Kxy3d[i3D] = Kxy_new * CUtoLUlength;
      Kxz3d[i3D] = Kxz_new * CUtoLUlength;
      Kyy3d[i3D] = Kyy_new * CUtoLUlength;
      Kyz3d[i3D] = Kyz_new * CUtoLUlength;
      Kzz3d[i3D] = Kzz_new * CUtoLUlength;
      Alpha3d[i3D] = Alpha_new;
      Betaux3d[i3D] = Betaux_new;
      Betauy3d[i3D] = Betauy_new;
      Betauz3d[i3D] = Betauz_new;
    }
  }

  // Some unit conversions and variable transformations 
  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++) {
        int i3D = i + nx*(j  + ny*k);

        P[i3D] = eps3d[i3D] * (((*Gamma) - 1.0) * rho3d[i3D]);

        // Apply floor 
        if (P[i3D]   < (*vacuum)*(*vacuum)) P[i3D] = (*vacuum)*(*vacuum);
        if (rho3d[i3D] < (*vacuum)) rho3d[i3D] = (*vacuum);

        
        theta[i3D] = 0;

      }
    }
  }

  free(x3d);
  free(y3d);
  free(z3d);

  magneticField(P, Sx, Sy, Sz, Bx, By, Bz, gxx3d, gxy3d, gxz3d, gyy3d, gyz3d, gzz3d, *tov_angle, *tov_Asize, *tov_facp, *vacuum, init_domain_x, init_domain_y, init_domain_z, init_patch_x, init_patch_y, init_patch_z, d_ghost_width, CUto10km, dx, nx, ny, nz);

  adm2bssn(eps3d, tau, rho3d, D, vx3d, 
           vy3d, vz3d, Sx, Sy, Sz, 
           P, Bx, By, Bz, gxx3d, 
           gxy3d, gxz3d, gyy3d, gyz3d, gzz3d,
           Kxx3d, Kxy3d, Kxz3d, Kyy3d, Kyz3d, 
           Kzz3d, Alpha3d, chi, trK, theta, 
           Gamh_x, Gamh_y, Gamh_z, *Gamma, Betaux3d, 
           Betauy3d, Betauz3d, dx, nx, ny, 
           nz);
}

void SAMRAI_External_Data::magneticField(double *P,   double *Sx,     double *Sy,     double *Sz, double *Bx,     double *By,     double *Bz, 
              double *gxx,   double *gxy,    double *gxz,
              double *gyy,   double *gyz,    double *gzz,
              double tov_angle, double tov_Asize, double tov_facp, double vacuum, double init_domain_x, double init_domain_y, double init_domain_z, int init_patch_x, int init_patch_y, int init_patch_z, int d_ghost_width, double CUto10KM, const double* dxA, int nx,        int ny,         int nz) {

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

void SAMRAI_External_Data::polint( double *xa, double *ya, int n, double x, double *y, double *dy )
{
    double *c = NULL;
    double *d = NULL;
    double den;
    double dif;
    double dift;
    double ho;
    double hp;
    double w;

    int i;
    int m;
    int ns;

    if( (c = (double *)malloc( n * sizeof( double ) )) == NULL ||
        (d = (double *)malloc( n * sizeof( double ) )) == NULL ) {
        fprintf( stderr, "polint error: allocating workspace\n" );
        fprintf( stderr, "polint error: setting y = 0 and dy = 1e9\n" );
        *y = 0.0;
        *dy = 1.e9;

        if( c != NULL )
            free( c );
        if( d != NULL )
            free( d );
        return;
    }

    ns = 0;
    dif = fabs(x-xa[0]);
    for( i = 0; i < n; ++i ) {
        dift = fabs( x-xa[i] );
        if( dift < dif ) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    *y = ya[ns];
    ns = ns-1;
    for( m = 0; m < n-1; ++m ) {
        for( i = 0; i < n-m-1; ++i ) {
            ho = xa[i]-x;
            hp = xa[i+m+1]-x;
            w = c[i+1]-d[i];
            den = ho-hp;
            if( den == 0 ) {
                fprintf( stderr, "polint error: den = 0\n" );
                fprintf( stderr, "polint error: setting y = 0 and dy = 1e9\n" );
                *y = 0.0;
                *dy = 1.e9;

                if( c != NULL )
                    free( c );
                if( d != NULL )
                    free( d );
                return;
            }
            den = w/den;
            d[i] = hp*den;
            c[i] = ho*den;
        }
        if( 2*(ns+1) < n-m-1 ) {
            *dy = c[ns+1];
        } else {
            *dy = d[ns];
            ns = ns-1;
        }
        *y = (*y)+(*dy);
    }


    if( c != NULL )
        free( c );
    if( d != NULL )
        free( d );
    return;
}
