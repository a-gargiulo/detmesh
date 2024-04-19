#define _USE_MATH_DEFINES
#include <math.h>

#include "diamond.h"

void arc_constraints(int n, double* x, double* fvec, int n_xtra_args, double* args) {
  double alpha = args[0] * M_PI / 180.0;
  double beta = args[1] * M_PI / 180.0;
  double l_d = args[2]; 
  double r = args[3]; 

  fvec[0] = (x[0] - x[2]) * (x[0] - x[2]) +
            (tan(alpha) * x[0] - x[3]) * (tan(alpha) * x[0] - x[3]) - r * r;
  fvec[1] =
      (x[1] - x[2]) * (x[1] - x[2]) +
      (tan(beta) * (l_d - x[1]) - x[3]) * (tan(beta) * (l_d - x[1]) - x[3]) -
      r * r;
  fvec[2] =
      x[0] * (x[2] - x[0]) + tan(alpha) * x[0] * (x[3] - tan(alpha) * x[0]);
  fvec[3] = (l_d - x[1]) * (x[1] - x[2]) -
            tan(beta) * (l_d - x[1]) * (tan(beta) * (l_d - x[1]) - x[3]);
}
