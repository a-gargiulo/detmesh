#include "diamond.h"

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
  if (f->size != 4) {
    fprintf(stderr, "Error: The size must be 4.\n");
    f->size = 0;
    return f;
  }
  f->data[0] = (x->data[0] - x->data[2]) * (x->data[0] - x->data[2]) +
         (leading_edge(x->data[0], geo->alpha) - x->data[3]) *
             (leading_edge(x->data[0], geo->alpha) - x->data[3]) -
         geo->r * geo->r;
  f->data[1] = (x->data[1] - x->data[2]) * (x->data[1] - x->data[2]) +
         (trailing_edge(x->data[1], geo->beta, geo->l_d) - x->data[3]) *
             (trailing_edge(x->data[1], geo->beta, geo->l_d) - x->data[3]) -
         geo->r * geo->r;
  f->data[2] = x->data[0] * (x->data[2] - x->data[0]) + leading_edge(x->data[0], geo->alpha) *
                                    (x->data[3] - leading_edge(x->data[0], geo->alpha));
  f->data[3] = (geo->l_d - x->data[1]) * (x->data[1] - x->data[2]) -
         trailing_edge(x->data[1], geo->beta, geo->l_d) *
             (trailing_edge(x->data[1], geo->beta, geo->l_d) - x->data[3]);
}


Vector* calculate_arc_parameters(const Vector* x_guess, const Geometry* geo) {

  Vector* x_old = create_vector(x_guess->size);
  Vector* x_new = create_vector(x_guess->size);
  Vector* rel_err = create_vector(x_guess->size);
  for (int i = 0; i < rel_err->size; ++i) {
    rel_err->data[i] = 1;
  }
  memcpy(x_new->data, x_guess->data, x_new->size*sizeof(double))

  while(vector_norm(rel_err)>1e-5){
    memcpy(x_old->data, x_new->data, x_old->size*sizeof(double)) 
    arc_constrain_functions(x_old,x_new, geo); 
    vector_subtract(x_new, x_old, rel_err);
    vector_division(rel_err, x_new, rel_err);
  }
  
  delete_vector(x_old);
  delete_vector(rel_err);

  return x_new
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
