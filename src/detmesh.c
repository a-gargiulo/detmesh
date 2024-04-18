#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <fsolve.h>
#include <math.h>
// #include <gmshc.h>

#include "diamond.h"
#include "parsing.h"
int main(int argc, char** argv) {
  Geometry geo;
  int n = 4;
  double* x = (double*)malloc(n * sizeof(double));
  double* fvec = (double*)calloc(n, sizeof(double));

  read_input(argv[1], &geo, x); 
  typedef void (*Fcn)(int, double*, double*);


  Fcn fcn = arc_constraints;
  double tol = 1e-5;
  fsolve(fcn, n, x, fvec, tol);

  for (int i = 0; i < n; ++i) {
    if (i == 3) {
      printf("%.12f\n", x[i]);
    } else {
      printf("%.12f, ", x[i]);
    }
  }
  free(x);
  free(fvec);
  return 0;
}
