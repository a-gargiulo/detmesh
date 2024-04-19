#ifndef DIAMOND_H
#define DIAMOND_H

typedef struct {
  double alpha;
  double beta;
  double l_d;
  double r;
} Geometry;

void arc_constraints(int n, double* x, double* fvec, int n_xtra_args,
                     double* args);

#endif  // DIAMOND_H
