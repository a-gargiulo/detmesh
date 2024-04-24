#include "diamond.h"

#define _USE_MATH_DEFINES
#define REL_TOL 1e-8

#include <fsolve.h>
#include <math.h>
#include <stdlib.h>

typedef void (*Fcn)(int, double*, double*, int, double*);

void arc_constraints(int n, double* x, double* fvec, int n_xtra_args,
                     double* args) {
  double alpha = args[0] * M_PI / 180.0;
  double beta = args[1] * M_PI / 180.0;
  double l = args[2];
  double r = args[3];

  fvec[0] = (x[0] - x[2]) * (x[0] - x[2]) +
            (tan(alpha) * x[0] - x[3]) * (tan(alpha) * x[0] - x[3]) - r * r;

  fvec[1] =
      (x[1] - x[2]) * (x[1] - x[2]) +
      (tan(beta) * (l - x[1]) - x[3]) * (tan(beta) * (l - x[1]) - x[3]) -
      r * r;

  fvec[2] =
      x[0] * (x[2] - x[0]) + tan(alpha) * x[0] * (x[3] - tan(alpha) * x[0]);

  fvec[3] = (l - x[1]) * (x[1] - x[2]) -
            tan(beta) * (l - x[1]) * (tan(beta) * (l - x[1]) - x[3]);
}

void calculate_arc_parameters(int n_arc_parameters, double* x_init,
                              Diamond* diamond) {
  Fcn fcn = arc_constraints;

  double rel_tol = REL_TOL;
  int n_xtra_args = 4;
  double* fvec = (double*)malloc(n_arc_parameters * sizeof(double));
  double* args = (double*)malloc(n_xtra_args * sizeof(double));
  args[0] = diamond->alpha;
  args[1] = diamond->beta;
  args[2] = diamond->l;
  args[3] = diamond->r;

  fsolve(fcn, n_arc_parameters, x_init, fvec, rel_tol, n_xtra_args, args);

  diamond->x1 = x_init[0];
  diamond->x2 = x_init[1];
  diamond->cx = x_init[2];
  diamond->cy = x_init[3];

  free(args);
  free(fvec);
}
