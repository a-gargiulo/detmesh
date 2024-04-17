#include "diamond.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utility.h"

double leading_edge(const double x, const float alpha) {
  return tan(alpha * M_PI / 180) * x;
}

double trailing_edge(const double x, const float beta, const double l_d) {
  return tan(beta) * (l_d - x);
}

void arc_constrain_functions(const double* x, const Geometry* geo, double* f) {
  f[0] = (x[0] - x[2]) * (x[0] - x[2]) +
         (leading_edge(x[0], geo->alpha) - x[3]) *
             (leading_edge(x[0], geo->alpha) - x[3]) -
         geo->r * geo->r;
  f[1] = (x[1] - x[2]) * (x[1] - x[2]) +
         (trailing_edge(x[1], geo->beta, geo->l_d) - x[3]) *
             (trailing_edge(x[1], geo->beta, geo->l_d) - x[3]) -
         geo->r * geo->r;
  f[2] = x[0] * (x[2] - x[0]) + leading_edge(x[0], geo->alpha) *
                                    (x[3] - leading_edge(x[0], geo->alpha));
  f[3] = (geo->l_d - x[1]) * (x[1] - x[2]) -
         trailing_edge(x[1], geo->beta, geo->l_d) *
             (trailing_edge(x[1], geo->beta, geo->l_d) - x[3]);
}

void calculate_arc_parameters(const double* x_guess, double* x, const Geometry* geo, int size) {
  double* x_old = (double*)malloc(size*sizeof(double));
  double* x_new = (double*)malloc(size*sizeof(double));
  while(){


  }
  free(x_old);
  free(x_new);
}


//   int MAXITER = 100;

//   double* x_old;
//   // memcpy(x_old, x_init,  size * sizeof(double));
//   double* x_new = x_init;
//   while (){

//     x_old = x_new;
//     x_new = arc_constrain_functions(x_old, f, geo);
//   }
// }
