#include "detmesh.h"

#include "diamond.h"
#include "mesh.h"
#include "parsing.h"

void run(int argc, char** argv) {
  Diamond diamond;
  double xInit[N_ARC_PARAMS]; 
  MeshConfig meshConfig;

  parse_input(argv[1], &diamond, xInit, &meshConfig);

  calculate_arc_parameters(xInit, &diamond);

  mesh_diamond(argc, argv, &diamond, &meshConfig);
}
