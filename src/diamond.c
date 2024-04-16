#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "geometry.h"

double leading_edge(double x, float alpha) {
  return tan(alpha * M_PI / 180) * x;
}

double trailing_edge(double x, float beta, double l_d) {
  return tan(beta) * (l_d - x);
}

void evaluate_constrain_equations(double* f, double* x, const Geometry* geom) {
  f[0] = pow(x[0] - x[2], 2) + pow(leading_edge(x[0], geom->alpha) - x[3], 2) -
         pow(geom->r, 2);
  f[1] = pow(x[1] - x[2], 2) +
         pow(trailing_edge(x[1], geom->beta, geom->l_d) - x[3], 2) -
         pow(geom->r, 2);
  f[2] = x[0] * (x[2] - x[0]) + leading_edge(x[0], geom->alpha) *
                                    (x[3] - leading_edge(x[0], geom->alpha));
  f[3] = (geom->l_d - x[1]) * (x[1] - x[2]) -
         trailing_edge(x[1], geom->beta, geom->l_d) *
             (trailing_edge(x[1], geom->beta, geom->l_d) - x[3])
}



Arc* calculate_arc(double* x, const Geometry* geom) {
  /* Arc* arc = (Arc*)malloc(sizeof(Arc)); */
  double f_old[4] = {0.0};

  while (f) {
    f_evaluate_constrain_equations(f, x, geom);
  }


}
