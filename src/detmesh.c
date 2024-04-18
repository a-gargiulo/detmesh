#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <fsolve.h>
#include <math.h>
// #include <gmshc.h>

// #include "diamond.h"
// #include "vector.h"
void arc(int n, double* x, double* fvec) {
  float alpha = 7.5 * M_PI / 180.0;
  float beta = 30 * M_PI / 180.0;
  double l_d = 0.7;
  double r = 0.0265;

  fvec[0] = (x[0] - x[2]) * (x[0] - x[2]) +
            (tan(alpha) * x[0] - x[3]) * (tan(alpha) * x[0] - x[3]) - r * r;
  fvec[1] =
      (x[1] - x[2]) * (x[1] - x[2]) +
      (tan(beta) * (l_d - x[1]) - x[3]) * (tan(beta) * (l_d - x[1]) - x[3]) -
      r * r;
  fvec[2] = x[0] * (x[2] - x[0]) + tan(alpha) * x[0] * (x[3] - tan(alpha) * x[0]);
  fvec[3] = (l_d - x[1]) * (x[1] - x[2]) -
         tan(beta) * (l_d - x[1]) * (tan(beta) * (l_d - x[1]) - x[3]);
}

int main() {
  typedef void (*Fcn)(int, double*, double*);

  double l_d = 0.7;
  int n = 4;
  double* x = (double*)malloc(n * sizeof(double));
  double* fvec = (double*)calloc(n, sizeof(double));
  x[0] = 0.5;
  x[1] = 0.55;
  x[2] = 0.52;
  x[3] = 0.01 * l_d;
  Fcn fcn = arc;
  double tol = 1e-5;
  fsolve(fcn, n, x, fvec, tol);

  for (int i = 0; i < n; ++i) {
    if (i == 3){
      printf("%.12f\n", x[i]);
    }
    else {
      printf("%.12f, ", x[i]);
    }
  }
  free(x);
  free(fvec);
  return 0;
}
