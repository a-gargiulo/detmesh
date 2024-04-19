#include "detmesh.h"
#include <stdio.h>
#include <stdlib.h>
#include "diamond.h"
#include "parsing.h"

void run(const char* input_file) {
  int n_arc_parameters = 4;
  Diamond* diamond = (Diamond*)malloc(sizeof(Diamond));
  double* x_init = (double*)malloc(n_arc_parameters * sizeof(double));

  read_input(input_file, diamond, x_init);
  calculate_arc_parameters(n_arc_parameters, x_init, diamond);

  free(x_init);
  free(diamond);

}
