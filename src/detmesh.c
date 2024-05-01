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
  size_t numEntityBlocks;
  size_t numEntityBlocksElem;

  Point* points = NULL;
  Curve* curves = NULL;
  Surface* surfaces = NULL;
  Node* nodes = NULL;
  Element* elements = NULL;

  int status;
  parse_input(argv[argc-1], &diamond, xInit, &meshConfig);

  calculate_arc_parameters(xInit, &diamond);

  status = mesh_diamond(argc, argv, &diamond, &meshConfig);
  // readGmsh("diamond.msh", &points, &numPoints, &curves, &numCurves, &surfaces, &numSurfaces, &nodes, &numEntityBlocks, &elements, &numEntityBlocksElem);

  ////CLEANUP
  //free(points);
  //free(curves);
  //free(surfaces);
  //for (size_t i = 0; i < numEntityBlocks; ++i) {
  //  free(nodes[i].nodeTags);
  //  free(nodes[i].x);
  //  free(nodes[i].y);
  //  free(nodes[i].z);
  //}
  //free(nodes);
  //for (size_t i = 0; i < numEntityBlocksElem; ++i) {
  //  free(elements[i].elementTags);
  //  free(elements[i].nodeTags);
  //}
  //free(elements);
}
