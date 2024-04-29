#include "detmesh.h"

#include <stdlib.h>

#include "diamond.h"
#include "mesh.h"
#include "parsing.h"
#include "converter.h"


void run(int argc, char** argv) {
  Diamond diamond;
  double xInit[N_ARC_PARAMS]; 
  MeshConfig meshConfig;
  Entities entities;

  parse_input(argv[argc-1], &diamond, xInit, &meshConfig);

  calculate_arc_parameters(xInit, &diamond);

  mesh_diamond(argc, argv, &diamond, &meshConfig);

  readGmsh("diamond.msh", &entities);

  free(entities.points);
  free(entities.curves);
  free(entities.surfaces);

}
