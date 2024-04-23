#include "mesh.h"

#define _USE_MATH_DEFINES

#include <gmshc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diamond.h"

#define deg2rad(deg) ((deg) * M_PI / 180.0)

int mesh_diamond(int argc, char** argv, Diamond* diamond) {
  // Mesh parameters
  double l_approach = 0.2;
  double l_wake = 0.2;
  double t_height = 0.1016;

  double l_blk1 = 0.08;
  double d_arc_u = 0.02;
  double d_arc_d = 0.02;
  double l_blk2 = 0.05;

  double lc = 5e-3;  // default mesh size
                     //
  // Local variables
  int ierr;

  // TO-DO: Include checks, e.g., l_blk_1 < l_approach
  // ...

  gmshInitialize(argc, argv, 1, 0, &ierr);

  if (gmshIsInitialized(&ierr)) {
    printf("SUCCESSFUL initialization of gmsh...\n");
  } else {
    printf("ERROR! Gmsh could not be initialize\n");
    return 1;
  }

  gmshModelAdd("detmesh", &ierr);

  // Fixed geometry points - 1 
  gmshModelGeoAddPoint(0, 0, 0, lc, 1, &ierr);           // 1
  gmshModelGeoAddPoint(l_approach, 0, 0, lc, 2, &ierr);  // 2
  gmshModelGeoAddPoint(l_approach + diamond->x1,
                       tan(deg2rad(diamond->alpha)) * diamond->x1, 0, lc / 5.0,
                       3,
                       &ierr);  // 3
  gmshModelGeoAddPoint(l_approach + diamond->cx, diamond->cy, 0, lc, 4,
                       &ierr);  // 4 Center
  gmshModelGeoAddPoint(
      l_approach + diamond->x2,
      tan(deg2rad(diamond->beta)) * (diamond->l_d - diamond->x2), 0, lc / 5.0,
      5, &ierr);                                                        // 5
  gmshModelGeoAddPoint(l_approach + diamond->l_d, 0, 0, lc, 6, &ierr);  // 6
  gmshModelGeoAddPoint(l_approach + diamond->l_d + l_wake, 0, 0, lc, 7,
                       &ierr);  // 7
  gmshModelGeoAddPoint(0, t_height, 0, lc, 8, &ierr);
  gmshModelGeoAddPoint(l_approach + diamond->l_d + l_wake, t_height, 0, lc, 9,
                       &ierr);


  // Block points - 100
  gmshModelGeoAddPoint(l_approach - l_blk1 / 2.0, 0, 0, lc, 100, &ierr);
  gmshModelGeoAddPoint(l_approach + l_blk1 / 2.0, tan(deg2rad(diamond->alpha)) * l_blk1 / 2.0, 0, lc, 101, &ierr);
  gmshModelGeoAddPoint(l_approach + diamond->x1 - d_arc_u, tan(deg2rad(diamond->alpha)) * (diamond->x1 - d_arc_u), 0, lc, 102, &ierr);
  gmshModelGeoAddPoint(l_approach + diamond->x2 + d_arc_d, tan(deg2rad(diamond->beta)) * (diamond->l_d - diamond->x2 - d_arc_d), 0, lc, 103, &ierr);
  gmshModelGeoAddPoint(l_approach + diamond->l_d - l_blk2 / 2.0, tan(deg2rad(diamond->beta)) * l_blk2/2.0, 0, lc, 104, &ierr);
  gmshModelGeoAddPoint(l_approach + diamond->l_d + l_blk2 / 2.0, 0, 0, lc, 105, &ierr);

  // gmshModelGeoAddPoint(l_approach - l_blk1 / 2.0, t_height, 0, lc, 106, &ierr);
  // gmshModelGeoAddPoint(l_approach + l_blk1 / 2.0, t_height, 0, lc, 107, &ierr);
  // gmshModelGeoAddPoint(l_approach + diamond->x1 - d_arc_u, t_height, 0, lc, 108, &ierr);
  // gmshModelGeoAddPoint(l_approach + diamond->x2 + d_arc_d, t_height, 0, lc, 109, &ierr);
  // gmshModelGeoAddPoint(l_approach + diamond->l_d - l_blk2 / 2.0, t_height, 0, lc, 110, &ierr);
  // gmshModelGeoAddPoint(l_approach + diamond->l_d + l_blk2 / 2.0, t_height, 0, lc, 111, &ierr);

  // Bottom Wall
  gmshModelGeoAddLine(1, 100, 1, &ierr);
  gmshModelGeoAddLine(100, 2, 2, &ierr);
  gmshModelGeoAddLine(2, 101, 3, &ierr);
  gmshModelGeoAddLine(101, 102, 4, &ierr);
  gmshModelGeoAddLine(102, 3, 5, &ierr);
  gmshModelGeoAddCircleArc(3, 4, 5, 6, 0, 0, 0, &ierr);

  gmshModelGeoAddLine(5, 103, 7, &ierr);
  gmshModelGeoAddLine(103, 104, 8, &ierr);
  gmshModelGeoAddLine(104, 6, 9, &ierr);
  gmshModelGeoAddLine(6, 105, 10, &ierr);
  gmshModelGeoAddLine(105, 7, 11, &ierr);


  // Top wall
  // gmshModelGeoAddLine(8, 9, 12, &ierr);
  // gmshModelGeoAddLine(8, 106, 11, &ierr);
  // gmshModelGeoAddLine(106, 107, 12, &ierr);
  // gmshModelGeoAddLine(107, 108, 13, &ierr);
  // gmshModelGeoAddLine(108, 109, 14, &ierr);
  // gmshModelGeoAddLine(109, 110, 15, &ierr);
  // gmshModelGeoAddLine(110, 111, 16, &ierr);
  // gmshModelGeoAddLine(111, 9, 17, &ierr);


  gmshModelGeoSynchronize(&ierr);

  double ratio = 1.05;
  double d_w = -1e-6;

  int dimTags[] = {1, 1, 1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9, 1, 10, 1, 11};
  size_t dimTags_n = 22;

  int* outDimTags = NULL;
  size_t outDimTags_n;

  size_t numElements_n = 100;
  int numElements[100];
  for (size_t i = 0; i < numElements_n; ++i) {
    numElements[i] = 1;
  }

  size_t heights_n = 100;
  double heights[100];
  for (size_t i = 0; i < heights_n; ++i) {
    if (i == 0) {
      heights[i] = d_w;
    } else {
      heights[i] = heights[i - 1] - (-d_w) * pow(ratio, i);
    }
  }

  int recombine = 1;
  int second = 0;
  int viewIndex = -1;

  gmshModelGeoExtrudeBoundaryLayer(
      dimTags, dimTags_n, &outDimTags, &outDimTags_n, numElements,
      numElements_n, heights, heights_n, recombine, second, viewIndex,
      &ierr);

  gmshFree(outDimTags);

  gmshModelGeoSynchronize(&ierr);
  gmshModelMeshGenerate(2, &ierr);
  // Figure out coordinates
  // int pts[] = {108, 109, 113, 117, 121, 125, 130, 134, 138, 142, 146, 150};
  int pts[] = {109, 117, 121, 134, 138, 146};
  size_t pts_n = 6;
  double* coord = NULL;
  size_t coord_n;
  for (size_t i = 0; i < pts_n; ++i) {
    gmshModelGetValue(0, pts[i], NULL, 0, &coord, &coord_n, &ierr);
    gmshModelGeoAddPoint(coord[coord_n-3], t_height, coord[coord_n-1], lc, 1000+i, &ierr);
    gmshFree(coord);
    coord = NULL;
  }

  gmshModelGeoAddLine(8, 1000, 100, &ierr);
  gmshModelGeoAddLine(1000, 1001, 101, &ierr);
  gmshModelGeoAddLine(1001, 1002, 102, &ierr);
  gmshModelGeoAddLine(1002, 1003, 103, &ierr);
  gmshModelGeoAddLine(1003, 1004, 104, &ierr);
  gmshModelGeoAddLine(1004, 1005, 105, &ierr);
  gmshModelGeoAddLine(1005, 9, 106, &ierr);
  // gmshModelGeoAddLine(1002, 1003, 14, &ierr);
  // gmshModelGeoAddLine(1003, 1004, 15, &ierr);

  double ratio2 = 1.1;
  double d_w2 = 1e-6;

  int dimTags2[] = {1, 100, 1, 101, 1, 102, 1, 103, 1, 104, 1, 105, 1, 106};
  size_t dimTags_n2 = 14;

  int* outDimTags2 = NULL;
  size_t outDimTags_n2;

  size_t numElements_n2 = 71;
  int numElements2[71];
  for (size_t i = 0; i < numElements_n2; ++i) {
    numElements2[i] = 1;
  }

  size_t heights_n2 = 71;
  double heights2[71];
  for (size_t i = 0; i < heights_n2; ++i) {
    if (i == 0) {
      heights2[i] = d_w2;
    } else {
      heights2[i] = heights2[i - 1] + d_w2 * pow(ratio2, i);
    }
  }

  int recombine2 = 1;
  int second2 = 0;
  int viewIndex2 = -1;

  gmshModelGeoExtrudeBoundaryLayer(dimTags2, dimTags_n2, &outDimTags2,
                                   &outDimTags_n2, numElements2,
                                   numElements_n2, heights2, heights_n2,
                                   recombine2, second2, viewIndex2, &ierr);
  gmshFree(outDimTags2);

  gmshModelGeoSynchronize(&ierr);
  gmshModelMeshGenerate(2, &ierr);

  gmshModelGeoAddLine(108, 1008, 200, &ierr);
  gmshModelGeoAddLine(109, 1009, 201, &ierr);
  gmshModelGeoAddLine(117, 1013, 202, &ierr);
  gmshModelGeoAddLine(121, 1017, 203, &ierr);
  gmshModelGeoAddLine(134, 1021, 204, &ierr);
  gmshModelGeoAddLine(138, 1025, 205, &ierr);
  gmshModelGeoAddLine(146, 1029, 206, &ierr);
  gmshModelGeoAddLine(150, 1033, 207, &ierr);


  gmshModelGeoSynchronize(&ierr);
  // Curve loops
  // gmshModelGeoAddLine(100, 1009, 888, &ierr);
  // gmshModelGeoAddLine(1, 1008, 889, &ierr);
  // int blk1[] = {888, 201, -889, -200};
  int blk1[] = {12, 201, -107, -200};
  // int blk1[] = {1, 201, -100, -200};
  // int blk1[] = {1, 888, -100, -889};
  int blk1_loops[] = {1};
  size_t blk1_n = 4; 
  size_t blk1_loops_n = 1;
  gmshModelGeoAddCurveLoop(blk1, blk1_n, 1, 0, &ierr);
  gmshModelGeoAddPlaneSurface(blk1_loops, blk1_loops_n, 300, &ierr);

  gmshModelGeoSynchronize(&ierr);
  gmshModelMeshGenerate(2, &ierr);

  // gmshModelGeoMeshSetTransfiniteCurve(888, 100, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(1, 100, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(200, 100, "Progression", 1, &ierr);
  // // gmshModelGeoMeshSetTransfiniteCurve(201, 100, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(889, 100, "Progression", 1, &ierr);
  // int cornerTags[] = {108, 109, 1009, 1008};
  // gmshModelGeoMeshSetTransfiniteSurface(1, "Left", cornerTags, 4, &ierr);
  // gmshModelGeoMeshSetRecombine(300, 1, 90, &ierr);
  // gmshModelGeoSynchronize(&ierr);
  // gmshModelMeshGenerate(2, &ierr);

  // gmshModelGeoAddLine(12, 33, 1000, &ierr);
  // gmshModelGeoSynchronize(&ierr);
  // int curveTags[] = {7, 11, 15, 19, 23, 999, -27, -1000};
  // size_t curveTags_n = 8;
  // gmshModelGeoAddCurveLoop(curveTags, curveTags_n, 1, 0, &ierr);
  // int wireTags[] = {1};
  // size_t wireTags_n = 1;
  // gmshModelGeoAddPlaneSurface(wireTags, wireTags_n, 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(1000, 100, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(999, 100, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(6, 200, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(1, 21, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(2, 61, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(3, 41, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(4, 41, "Progression", 1, &ierr);
  // gmshModelGeoMeshSetTransfiniteCurve(5, 40, "Progression", 1, &ierr);
  // // gmshModelGeoMeshSetRecombine(2, 1, 90, &ierr);
  // int cornerTags[] = {12, 30, 34, 33};
  // gmshModelGeoMeshSetTransfiniteSurface(1, "Left", &cornerTags, 4, &ierr);
  // gmshModelGeoMeshSetRecombine(2, 1, 90, &ierr);

  // int loopTags[] = {7, 11, 15, 19, 23, 101, -27, -100};
  // gmshModelGeoAddCurveLoop(loopTags, 8, 1, 0 ,&ierr);
  // int wireTags[] = {1};
  // gmshModelGeoAddPlaneSurface(wireTags, 1, 1, &ierr);
  gmshWrite("diamond.msh", &ierr);

  int gui = 1;
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-nopopup")) {
      gui = 0;
      break;
    }
  }
  if (gui) gmshFltkRun(&ierr);

  gmshFinalize(&ierr);
  return 0;
  // free(outDimTags);
}
