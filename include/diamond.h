#ifndef DIAMOND_H
#define DIAMOND_H

#define REL_TOL 1e-8

typedef struct {
  float alpha;
  float beta;
  double l_d;
  double r;
} Geometry;

void arc_constraints(int n, double* x, double* fvec);

#endif  // DIAMOND_H
