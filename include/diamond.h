#ifndef DIAMOND_H
#define DIAMOND_H

typedef struct {
  float alpha;
  float beta;
  double l_d;
  double r;
} Geometry;

double leading_edge(double x, float alpha);

double trailing_edge(double x, float beta, double l_d);

void arc_constrain_functions(double* x, double* f, const Geometry* geo);

#endif  // DIAMOND_H
