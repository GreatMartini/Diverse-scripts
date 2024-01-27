#include "Functions.h"
#include <string>


#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define SIGN(X) (((X) > 0) - ((X) < 0))

using namespace external;

#define vector(v, i, j, k) (v)[i+ilast*(j+jlast*(k))]


