#ifndef DIAMOND_H
#define DIAMOND_H

#define N_ARC_PARAMS 4

typedef struct {
  double alpha;
  double beta;
  double l;
  double r;
  double x1;
  double x2;
  double cx;
  double cy;
} Diamond;

void arc_constraints(int n, double* x, double* fvec, int nXtraArgs,
                     double* xtraArgs);

void calculate_arc_parameters(double* xInit, Diamond* diamond);

#endif  // DIAMOND_H
