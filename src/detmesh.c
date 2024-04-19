#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <fsolve.h>
#include <math.h>
// #include <gmshc.h>

#include "diamond.h"
#include "parsing.h"

int main(int argc, char** argv) {
  int n = 4;
  Geometry geo;
  double* x = (double*)malloc(n * sizeof(double));
  read_input(argv[1], &geo, x);
  double* fvec = (double*)calloc(n, sizeof(double));

  typedef void (*Fcn)(int, double*, double*, int, double*);

  Fcn fcn = arc_constraints;
  double tol = 1e-5;
  int n_xtra_args = 4;
  double* args = (double*)malloc(n_xtra_args * sizeof(double)); 
  args[0] = geo.alpha;
  args[1] = geo.beta;
  args[2] = geo.l_d;
  args[3] = geo.r;
  fsolve(fcn, n, x, fvec, tol, n_xtra_args, args);

  for (int i = 0; i < n; ++i) {
    if (i == 3) {
      printf("%.12f\n", x[i]);
    } else {
      printf("%.12f, ", x[i]);
    }
  }
  free(x);
  free(fvec);
  free(args);
  return 0;
}
