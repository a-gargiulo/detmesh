#include "diamond.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vector.h"

double leading_edge(const double x, const float alpha) {
  return tan(alpha * M_PI / 180) * x;
}

double trailing_edge(const double x, const float beta, const double l_d) {
  return tan(beta) * (l_d - x);
}

void arc_constrain_functions(const Vector* x, Vector* f, const Geometry* geo) {
  // if (f->size != 4) {
  //   fprintf(stderr, "ERROR! The equation vector must be of size 4!\n");
  //   return -1;
  // }
  f->data[0] = (x->data[0] - x->data[2]) * (x->data[0] - x->data[2]) +
               (leading_edge(x->data[0], geo->alpha) - x->data[3]) *
                   (leading_edge(x->data[0], geo->alpha) - x->data[3]) -
               geo->r * geo->r;
  f->data[1] =
      (x->data[1] - x->data[2]) * (x->data[1] - x->data[2]) +
      (trailing_edge(x->data[1], geo->beta, geo->l_d) - x->data[3]) *
          (trailing_edge(x->data[1], geo->beta, geo->l_d) - x->data[3]) -
      geo->r * geo->r;
  f->data[2] = x->data[0] * (x->data[2] - x->data[0]) +
               leading_edge(x->data[0], geo->alpha) *
                   (x->data[3] - leading_edge(x->data[0], geo->alpha));
  f->data[3] =
      (geo->l_d - x->data[1]) * (x->data[1] - x->data[2]) -
      trailing_edge(x->data[1], geo->beta, geo->l_d) *
          (trailing_edge(x->data[1], geo->beta, geo->l_d) - x->data[3]);
}

Vector* calculate_arc_parameters(const Vector* x_guess, const Geometry* geo) {
  Vector* x_old = construct_vector(x_guess->size);
  Vector* x_new = construct_vector(x_guess->size);
  Vector* rel_err = construct_vector(x_guess->size);

  for (int i = 0; i < rel_err->size; ++i) {
    rel_err->data[i] = 1;
  }
  memcpy(x_new->data, x_guess->data, x_new->size * sizeof(double));

  int counter = 0;
  while (counter < MAX_ITER) {
    memcpy(x_old->data, x_new->data, x_old->size * sizeof(double));
    for (int i = 0; i < 4; ++i) {
      printf("%1.4f, ", x_new->data[i]);
    }
    printf("\n");
    arc_constrain_functions(x_old, x_new, geo);
    vector_subtract(x_new, x_old, rel_err);
    vector_divide(rel_err, x_new, rel_err);
    if (vector_norm(rel_err) < REL_TOL) {
      printf("CONVERGED!\n");
      break;
    }
    counter++;
    if (counter >= MAX_ITER){
      printf("MAXIMUM ITER REACHED!\n");
    }
  }

  delete_vector(x_old);
  delete_vector(rel_err);

  return x_new;
}
