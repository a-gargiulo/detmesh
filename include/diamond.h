#ifndef DIAMOND_H
#define DIAMOND_H

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

void arc_constraints(int n, double* x, double* fvec, int n_xtra_args,
                     double* args);


void calculate_arc_parameters(int n_arc_parameters, double* x_init, Diamond* diamond);


#endif  // DIAMOND_H
