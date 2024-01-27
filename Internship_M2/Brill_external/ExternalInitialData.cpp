// En este codigo gurardamos partes de los ejemplos que no se utilizaron, peros se guardaron como forma de referencia
#include <iostream>
#include <cmath>
#include <cstdarg>
#include <sys/stat.h>
#include <cstdlib>
#include "hdf5.h"

#include "ExternalInitialData.h"
#include "Functions.h"

using namespace std;
// Funcion que lee las variables
void SAMRAI_External_Data::external_loadData(int vars, ...) {
//Cargamos las variables del problema, Solo nos interesa Phi y las variables computacionales. No se supo como operar bien con
// las variables computacionales ya que el codigo no funcionaba bien al usarlas.
  va_list args; //Inciamos la lista de variables
  va_start(args, vars); 
  double* Phi = va_arg(args, double *); // Phi: el logaritmo del factor conforme
  double* dx = va_arg(args, double *); // El paso en x
  // Variables de computacionales de SAMRAI
  double init_domain_x = va_arg(args, double); 
  double init_domain_y = va_arg(args, double);
  double init_domain_z = va_arg(args, double);
  int init_patch_x = va_arg(args, int);
  int init_patch_y = va_arg(args, int);
  int init_patch_z = va_arg(args, int);
  int d_ghost_width = va_arg(args, int); // Tamano de la zona ghost
  int nx = va_arg(args, int); // Numero de celdas en x
  int ny = va_arg(args, int); // Numero de celdas en y
  int nz = va_arg(args, int); // Numero de celdas en z
  va_end(args); // Terminamos la lista de variables

  //Get position and data from initialization files
  std::string fileName; // Declaramos la variable nombre del archivo
  // Se declaran parametros intrinsecos a los archivos hdf5
  hid_t       file, space, dset;
  herr_t      status;

  int dataspace, rank; //Variables introducidas inicialmente para guardar el tamano del archivo

  fileName = "phi.h5"; // Nombre del archivp con datos a cargar
  file = H5Fopen (fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); // Abrimos el archivo
  dset = H5Dopen (file, "Phi", H5P_DEFAULT); // Leemos los datos de Phi


  double* data_phi = new double[100*100*100]; // Creamos un arreglo unidimensional

  //Optenemos los datos
  dset = H5Dopen (file, "Phi", H5P_DEFAULT);
  status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_phi);
  status = H5Dclose (dset);

  //Cerramos el archivo
  status = H5Fclose (file);

  // Numero total de datos
  int nd = nx * ny * nz;

  // rescale the coordinates
  for (int k=0; k<nx; k++) { // Loop sobre x
    for (int j=0; j<ny; j++) { // Loop sobre y
      for (int i=0; i<nz; i++) { // Loop sobre z
        double x3d; // Variable no usada pero declarada en los ejemplos
        double y3d; // Variable no usada pero declarada en los ejemplos
        double z3d; // Variable no usada pero declarada en los ejemplos
        x3d = (init_domain_x + ((i + init_patch_x) - d_ghost_width) * dx[0]); // Variable no usada pero declarada en los ejemplos
        y3d = (init_domain_y + ((j + init_patch_y) - d_ghost_width) * dx[1]); // Variable no usada pero declarada en los ejemplos
        z3d = (init_domain_z + ((k + init_patch_z) - d_ghost_width) * dx[2]); // Variable no usada pero declarada en los ejemplos
        int i3D = i + nx*(j  + ny*k); //Calculamos el indice tridimencional (Ver los comentarios al final del codigo para mas informacion)
        if(i3D<34671){ // SAMRAI se salta 34671 antes de empezar a leer los datos
          Phi[i3D]=0; // Introducimos Phi=0 en esta region
          }
        else{
          if(i<=102 && j<=102 && k<=102){ // SAMRAI se salta varios datos cuando se rebasa de los 100 datos (+ 3 puntos de ghost zone)
            Phi[i3D] = data_phi[i3D-34671-7*(j-3)-1449*(k-3)]; // Asignamos los valores de Phi
          } // SAMRAI genera un dominio de n+1 datos en cada eje entonces se trataron de copiar los datos del punto n al punto n+1 para completar el dominio artificialmenete
          else if(i>102 && j<=102 && k<=102){ // Esta parte es para extender artificialmente los datos de Phi al dominio completo de n+1
            Phi[i3D] = data_phi[i3D-34671-7*(j-3)-1449*(k-3)-1];
          }
          else if(i<=102 && j>102 && k<=102){// Esta parte es para extender artificialmente los datos de Phi al dominio completo de n+1
            Phi[i3D] = data_phi[i3D-34671-7*(j-4)-1449*(k-3)];
          }
          else if(i>102 && j>102 && k<=102){// Esta parte es para extender artificialmente los datos de Phi al dominio completo de n+1
            Phi[i3D] = data_phi[i3D-34671-7*(j-4)-1449*(k-3)-1];
          }
          else if(k>102){
            Phi[i3D] = data_phi[i3D-34671-7*(j-3)-1449*(k-4)];//data_phi[i3D-34671-6*(j-3)-1248*(k-3)];
          }
          else{
            Phi[i3D] = 0;
          }
        }
      }
    }
  }

  free(data_phi); // Liberamos los datos
}
// SAMRAI va a leer los datos de la siguiente manera:
/*
- Los primeros datos que lee el indice tridimensional son sobre x para los primeros puntos de y y de z,
empezando por la parte inferior del dominio. Luego pasa a llenar los datos sobre x para el segundo punto en y
y el primero de z etc... 

- SAMRAI tiene una zona fantasma de 3 puntos sobre cada eje.

- La lectura de datos se salta inicialmente 34671 datos sobre el eje z. Al pasar a un nuevo punto en z se salta
1449 datos y al pasar a un nuevo punto de y se salta 7 datos. No se entiende bien por que el programa lee los datos 
de esta manera

- Para probar la lectura de datos recomendamos generar un arreglo tridimensional de datos del 0 a n y, correr el programa
y analizar como se inicializan los datos. Este arreglo se puede generar en el codigo python.
*/

