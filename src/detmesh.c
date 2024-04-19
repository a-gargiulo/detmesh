#include "detmesh.h"

#include <stdlib.h>

#include "diamond.h"
#include "mesh.h"
#include "parsing.h"

void run(int argc, char** argv) {
  int n_arc_parameters = 4;
  Diamond* diamond = (Diamond*)malloc(sizeof(Diamond));
  double* x_init = (double*)malloc(n_arc_parameters * sizeof(double));

  read_input(argv[1], diamond, x_init);
  calculate_arc_parameters(n_arc_parameters, x_init, diamond);
  mesh_diamond(argc, argv, diamond);

  free(x_init);
  free(diamond);
}
