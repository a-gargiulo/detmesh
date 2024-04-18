#define _USE_MATH_DEFINES
#include <math.h>

#include "diamond.h"

void arc_constraints(int n, double* x, double* fvec) {
  float alpha = 7.5 * M_PI / 180.0;
  float beta = 30 * M_PI / 180.0;
  double l_d = 0.7;
  double r = 0.0265;

  // float alpha = geo->alpha;
  // float beta = geo->beta;
  // double l_d = geo->l_d;
  // double r = geo->r;

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
