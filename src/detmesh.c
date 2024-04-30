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

  // for (size_t i = 0; i < numPoints; ++i) {
  //   printf("%d %.12lf %.12lf %.12lf\n", points[i].tag, points[i].x, points[i].y, points[i].z);
  // }

  // for (size_t i = 0; i < numCurves; ++i) {
  //   printf("%d %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf 0 2 %d %d\n", curves[i].tag, curves[i].minX, curves[i].minY, curves[i].minZ, curves[i].maxX, curves[i].maxY, curves[i].maxZ, curves[i].tagsBoundingPoints[0], curves[i].tagsBoundingPoints[1]);
  // }

  // for (size_t i = 0; i < numSurfaces; ++i) {
  //   printf("%d %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf 0 4 %d %d %d %d\n", surfaces[i].tag, surfaces[i].minX, surfaces[i].minY, surfaces[i].minZ, surfaces[i].maxX, surfaces[i].maxY, surfaces[i].maxZ, surfaces[i].tagsBoundingCurves[0], surfaces[i].tagsBoundingCurves[1],surfaces[i].tagsBoundingCurves[2], surfaces[i].tagsBoundingCurves[3]);
  // }

  free(points);
  free(curves);
  free(surfaces);
  free(nodes);

}
