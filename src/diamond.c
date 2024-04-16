#include "diamond.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double leading_edge(double x, float alpha) {
  return tan(alpha * M_PI / 180) * x;
}

double trailing_edge(double x, float beta, double l_d) {
  return tan(beta) * (l_d - x);
}

void arc_constrain_functions(double* x, double* f, const Geometry* geo) {
  f[0] = pow(x[0] - x[2], 2) + pow(leading_edge(x[0], geo->alpha) - x[3], 2) -
         pow(geo->r, 2);
  f[1] = pow(x[1] - x[2], 2) +
         pow(trailing_edge(x[1], geo->beta, geo->l_d) - x[3], 2) -
         pow(geo->r, 2);
  f[2] = x[0] * (x[2] - x[0]) + leading_edge(x[0], geo->alpha) *
                                    (x[3] - leading_edge(x[0], geo->alpha));
  f[3] = (geo->l_d - x[1]) * (x[1] - x[2]) -
         trailing_edge(x[1], geo->beta, geo->l_d) *
             (trailing_edge(x[1], geo->beta, geo->l_d) - x[3]);
}
