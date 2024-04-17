#ifndef DIAMOND_H
#define DIAMOND_H

#define MAXITER 100


typedef struct {
  float alpha;
  float beta;
  double l_d;
  double r;
} Geometry;

double leading_edge(const double x, const float alpha);

double trailing_edge(const double x, const float beta, const double l_d);

void arc_constrain_functions(const double* x, const Geometry* geo, double* f);

#endif  // DIAMOND_H
