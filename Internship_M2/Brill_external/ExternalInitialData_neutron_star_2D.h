
using namespace std;

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


};



