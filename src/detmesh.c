#include <stdio.h>
#include <gmshc.h>
#include "diamond.h" 
#include "vector.h"

int main(int argc, char **argv)
{
  int vec_size = 4;
  Geometry geo; 
  geo.alpha = 7.5; 
  geo.beta = 30;
  geo.l_d = 0.7;
  geo.r = 0.0265;

  Vector* x_guess = construct_vector(vec_size); 
  x_guess->data[0] = 0.4*0.7;
  x_guess->data[1] = 0.9*0.7;
  x_guess->data[2] = 0.5*0.7;
  x_guess->data[3] = 0.01*0.7;

  Vector* arc_parameters = calculate_arc_parameters(x_guess, &geo);

  delete_vector(arc_parameters);
  return 0;
}
