#ifndef DIAMOND_H
#define DIAMOND_H

#include "vector.h"

#define MAX_ITER 100000
#define REL_TOL 1e-8

typedef struct {
  float alpha;
  float beta;
  double l_d;
  double r;
} Geometry;

double leading_edge(const double x, const float alpha);

double trailing_edge(const double x, const float beta, const double l_d);

void arc_constrain_functions(const Vector* x, Vector* f, const Geometry* geo);

Vector* calculate_arc_parameters(const Vector* x_guess, const Geometry* geo);

#endif  // DIAMOND_H
