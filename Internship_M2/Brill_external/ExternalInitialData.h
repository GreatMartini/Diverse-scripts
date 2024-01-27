// Cabecera de la funcion de lectura de datos
using namespace std;

class SAMRAI_External_Data
{

public:

   static void external_loadData(int vars, ...);

private:

  static void adm2bssn(double *phi, const double* dxA, int nx, int ny, int nz);


};



