#include "detmesh.h"

#include <stdio.h>
#include <stdlib.h>

#include "diamond.h"
#include "mesh.h"
#include "parsing.h"
#include "converter.h"


void run(int argc, char** argv) {
  Diamond diamond;
  double xInit[N_ARC_PARAMS]; 
  MeshConfig meshConfig;

  size_t numPoints, numCurves, numSurfaces;

  Point* points = NULL;
  Curve* curves = NULL;
  Surface* surfaces = NULL;
  Node* nodes = NULL;

  int status;

  parse_input(argv[argc-1], &diamond, xInit, &meshConfig);

  calculate_arc_parameters(xInit, &diamond);

  status = mesh_diamond(argc, argv, &diamond, &meshConfig);
  readGmsh("diamond.msh", &points, &numPoints, &curves, &numCurves, &surfaces, &numSurfaces, &nodes);

  for (size_t i = 0; i < numPoints; ++i) {
    printf("%d %.12lf %.12lf %.12lf\n", points[i].tag, points[i].x, points[i].y, points[i].z);
  }
  free(points);
  free(curves);
  free(surfaces);
  free(nodes);

}
