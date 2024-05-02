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

  int numPoints, numCurves, numSurfaces;
  int numEntityBlocks;
  int numEntityBlocksElem;
  int numNodes;

  Point* points = NULL;
  Curve* curves = NULL;
  Surface* surfaces = NULL;
  Node* nodes = NULL;
  Element* elements = NULL;

  int status;
  parse_input(argv[argc-1], &diamond, xInit, &meshConfig);

  calculate_arc_parameters(xInit, &diamond);

  status = mesh_diamond(argc, argv, &diamond, &meshConfig);
  readGmsh("diamond.msh", &points, &numPoints, &curves, &numCurves, &surfaces, &numSurfaces, &nodes, &numEntityBlocks, &elements, &numEntityBlocksElem, &numNodes);

  writeFluent("diamond_fluent.msh", nodes, numEntityBlocks, numNodes, &diamond, &meshConfig);

  ////CLEANUP
  free(points);
  free(curves);
  free(surfaces);
  for (int i = 0; i < numEntityBlocks; ++i) {
    free(nodes[i].nodeTags);
    free(nodes[i].x);
    free(nodes[i].y);
    free(nodes[i].z);
  }
  free(nodes);
  for (int i = 0; i < numEntityBlocksElem; ++i) {
    free(elements[i].elementTags);
    free(elements[i].nodeTags);
  }
  free(elements);
}
