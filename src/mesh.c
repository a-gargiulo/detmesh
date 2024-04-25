#include "mesh.h"

#define _USE_MATH_DEFINES
#include <gmshc.h>
#include <math.h>
#include <stdio.h> 
#include <string.h>

#include "diamond.h"

#define DEFAULT_CELL_SIZE 5e-3
#define deg2rad(deg) ((deg) * M_PI / 180.0)


int mesh_diamond(int argc, char** argv, Diamond* diamond, MeshConfig* meshConfig) {
  double dAlpha = deg2rad(diamond->alpha);
  double dBeta = deg2rad(diamond->beta);
  double dL = diamond->l;
  double dX1 = diamond->x1;
  double dX2 = diamond->x2;
  double dCx = diamond->cx;
  double dCy = diamond->cy;

  double sUp = meshConfig->sUp; 
  double sDown = meshConfig->sDown; 
  double sBlkUp = meshConfig->sBlkUp;
  double sBlkDown = meshConfig->sBlkDown;
  double sArcUp = meshConfig->sArcUp;
  double sArcDown = meshConfig->sArcDown;
  double tHeight = meshConfig->tHeight;

  int ierr;

  // TO-DO: Include checks, e.g., l_blk_1 < sUp

  gmshInitialize(argc, argv, 1, 0, &ierr);

  if (gmshIsInitialized(&ierr)) {
    printf("SUCCESSFUL initialization of gmsh...\n");
  } else {
    printf("ERROR! Gmsh could not be initialize\n");
    return 1;
  }

  gmshModelAdd("detmesh", &ierr);


  // GEOMETRY
  // -------- 
  double sPts[] = {
    0, 0,
    sUp - sBlkUp / 2.0, 0, 
    sUp, 0,
    sUp + sBlkUp / 2.0, tan(dAlpha) * sBlkUp / 2.0, 
    sUp + dX1 - sArcUp,  tan(dAlpha) * (dX1 - sArcUp),
    sUp + dX1, tan(dAlpha) * dX1,
    sUp + dCx, dCy,
    sUp + dX2, tan(dBeta) * (dL - dX2),
    sUp + dX2 + sArcDown,  tan(dBeta) * (dL - dX2 - sArcDown),
    sUp + dL - sBlkDown / 2.0, tan(dBeta) * sBlkDown / 2.0,
    sUp + dL, 0,
    sUp + dL + sBlkDown / 2.0, 0,
    sUp + dL + sDown, 0
  };

  size_t sPts_n = 26; 

  // Add bottom points (1 - 13)
  int tagCount = 0;
  for (size_t i = 0; i < sPts_n; i += 2) {
    gmshModelGeoAddPoint(sPts[i], sPts[i + 1], 0, DEFAULT_CELL_SIZE, ++tagCount, &ierr);
  }

  // Add bottom lines (1 - 11) 
  tagCount = 0;
  for (size_t i = 1; i < sPts_n / 2; ++i) {
    if (i == 6)
      gmshModelGeoAddCircleArc(i, i + 1, i + 2, ++tagCount, 0, 0, 0, &ierr);
    else if (i == 7)
      continue;
    else
      gmshModelGeoAddLine(i, i + 1, ++tagCount, &ierr);
  }

  // Extrude bottom boundary layer
  double sBlRatio = 1.05;
  double sBlW = -1e-6;

  int sBlTags[] = {1, 1, 1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9, 1, 10, 1, 11};
  size_t sBlTags_n = 22;

  int* sBlOutTags = NULL;
  size_t sBlOutTags_n;

  size_t sBlSteps = 100;
  int sBlElementsPerStep[100];
  for (size_t i = 0; i < sBlSteps; ++i) {
    sBlElementsPerStep[i] = 1;
  }

  size_t sBlHeights_n = 100;
  double sBlHeights[100];
  for (size_t i = 0; i < sBlHeights_n; ++i) {
    if (i == 0) {
      sBlHeights[i] = sBlW;
    } else {
      sBlHeights[i] = sBlHeights[i - 1] - (-sBlW) * pow(sBlRatio, i);
    }
  }

  int recombine = 1;
  int second = 0;
  int viewIndex = -1;

  // c12 p(16 17), c16 p(17 21), c20 p(21 25), c24 p(25 29), c28 p(29, 33), c32 p(33 38), c36 p(38 42), c40 p(42 46), c44 p(46 50), c48 p(50 54), c52 p(54 58)
  gmshModelGeoExtrudeBoundaryLayer(
      sBlTags, sBlTags_n, &sBlOutTags, &sBlOutTags_n, sBlElementsPerStep,
      sBlSteps, sBlHeights, sBlHeights_n, recombine, second, viewIndex,
      &ierr);

  gmshModelGeoSynchronize(&ierr);
  gmshModelMeshGenerate(2, &ierr);
 
  // index family 100
  int connectingBlPts[] = {17, 25, 29, 42, 46, 54};
  size_t connectingBlPts_n = 6;
  double* coord = NULL;
  size_t coord_n;
  for (size_t i = 0; i < connectingBlPts_n; ++i) {
    gmshModelGetValue(0, connectingBlPts[i], NULL, 0, &coord, &coord_n, &ierr);
    gmshModelGeoAddPoint(coord[coord_n-3], tHeight, coord[coord_n-1], DEFAULT_CELL_SIZE, 101 + i, &ierr);
    gmshFree(coord);
    coord = NULL;
  }
  gmshModelGeoAddPoint(0, tHeight, 0, DEFAULT_CELL_SIZE, 100, &ierr);
  gmshModelGeoAddPoint(sUp + dL + sDown, tHeight, 0, DEFAULT_CELL_SIZE, 107, &ierr);

  gmshModelMeshClear(0, 0, &ierr);
  gmshModelGeoSynchronize(&ierr);

//   gmshModelGeoMeshSetTransfiniteCurve(1, 100, "Progression", 1, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(2, 50, "Progression", 1, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(3, 51, "Progression", 1, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(4, 100, "Progression", 1, &ierr);

//   // Top wall
//   // gmshModelGeoAddLine(8, 9, 12, &ierr);
//   // gmshModelGeoAddLine(8, 106, 11, &ierr);
//   // gmshModelGeoAddLine(106, 107, 12, &ierr);
//   // gmshModelGeoAddLine(107, 108, 13, &ierr);
//   // gmshModelGeoAddLine(108, 109, 14, &ierr);
//   // gmshModelGeoAddLine(109, 110, 15, &ierr);
//   // gmshModelGeoAddLine(110, 111, 16, &ierr);
//   // gmshModelGeoAddLine(111, 9, 17, &ierr);


  // gmshModelGeoSynchronize(&ierr);


//   gmshModelGeoSynchronize(&ierr);
//   gmshModelMeshGenerate(2, &ierr);

//   gmshModelGeoAddLine(8, 1000, 100, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(100, 100, "Progression", 1, &ierr);
//   gmshModelGeoAddLine(1000, 1001, 101, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(101, 100, "Progression", 1, &ierr);
//   gmshModelGeoAddLine(1001, 1002, 102, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(102, 100, "Progression", 1, &ierr);
//   gmshModelGeoAddLine(1002, 1003, 103, &ierr);
//   gmshModelGeoAddLine(1003, 1004, 104, &ierr);
//   gmshModelGeoAddLine(1004, 1005, 105, &ierr);
//   gmshModelGeoAddLine(1005, 9, 106, &ierr);
//   // gmshModelGeoAddLine(1002, 1003, 14, &ierr);
//   // gmshModelGeoAddLine(1003, 1004, 15, &ierr);

//   double ratio2 = 1.1;
//   double d_w2 = 1e-6;

//   int dimTags2[] = {1, 100, 1, 101, 1, 102, 1, 103, 1, 104, 1, 105, 1, 106};
//   size_t dimTags_n2 = 14;

//   int* outDimTags2 = NULL;
//   size_t outDimTags_n2;

//   size_t numElements_n2 = 71;
//   int numElements2[71];
//   for (size_t i = 0; i < numElements_n2; ++i) {
//     numElements2[i] = 1;
//   }

//   size_t heights_n2 = 71;
//   double heights2[71];
//   for (size_t i = 0; i < heights_n2; ++i) {
//     if (i == 0) {
//       heights2[i] = d_w2;
//     } else {
//       heights2[i] = heights2[i - 1] + d_w2 * pow(ratio2, i);
//     }
//   }

//   int recombine2 = 1;
//   int second2 = 0;
//   int viewIndex2 = -1;

//   gmshModelGeoExtrudeBoundaryLayer(dimTags2, dimTags_n2, &outDimTags2,
//                                    &outDimTags_n2, numElements2,
//                                    numElements_n2, heights2, heights_n2,
//                                    recombine2, second2, viewIndex2, &ierr);
//   gmshFree(outDimTags2);

//   gmshModelGeoSynchronize(&ierr);

//   gmshModelGeoAddLine(108, 1008, 200, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(200, 50, "Progression", 1, &ierr);
//   gmshModelGeoAddLine(109, 1009, 201, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(201, 50, "Progression", 1, &ierr);
//   gmshModelGeoAddLine(117, 1013, 202, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(202, 50, "Progression", 1, &ierr);
//   gmshModelGeoAddLine(121, 1017, 203, &ierr);
//   gmshModelGeoMeshSetTransfiniteCurve(203, 50, "Progression", 1, &ierr);
//   gmshModelGeoAddLine(134, 1021, 204, &ierr);
//   gmshModelGeoAddLine(138, 1025, 205, &ierr);
//   gmshModelGeoAddLine(146, 1029, 206, &ierr);
//   gmshModelGeoAddLine(150, 1033, 207, &ierr);


//   // Curve loops
//   // gmshModelGeoAddLine(108, 109, 888, &ierr);
//   // gmshModelGeoAddLine(1008, 1009, 889, &ierr);
//   // int blk1[] = {888, 201, -889, -200};
//   int blk1[] = {12, 201, -107, -200};
//   // int blk1[] = {1, 201, -100, -200};
//   // int blk1[] = {1, 888, -100, -889};
//   int blk1_loops[] = {1};
//   size_t blk1_n = 4; 
//   size_t blk1_loops_n = 1;
//   gmshModelGeoAddCurveLoop(blk1, blk1_n, 1, 0, &ierr);
//   gmshModelGeoAddPlaneSurface(blk1_loops, blk1_loops_n, 300, &ierr);


//   int blk2[] = {16, 20, 202,-111, -201};
//   size_t blk2_n = 5;
//   int blk2_loops[]={2};
//   size_t blk2_loops_n = 1;
//   gmshModelGeoAddCurveLoop(blk2, blk2_n, 2, 0, &ierr);
//   gmshModelGeoAddPlaneSurface(blk2_loops, blk2_loops_n, 301, &ierr);

//    int blk3[] = {24, 203,-115, -202};
//   size_t blk3_n = 4;
//   int blk3_loops[]={3};
//   size_t blk3_loops_n = 1;
//   gmshModelGeoAddCurveLoop(blk3, blk3_n, 3, 0, &ierr);
//   gmshModelGeoAddPlaneSurface(blk3_loops, blk3_loops_n, 302, &ierr);

//   gmshModelGeoSynchronize(&ierr);
//  int cornerTags2[] = {109, 117, 1013, 1009};
//   gmshModelGeoMeshSetTransfiniteSurface(301, "Left", cornerTags2, 4, &ierr);
//   gmshModelGeoMeshSetRecombine(2,301, 90, &ierr);


//   int cornerTags3[] = {117,121,1017, 1013};
//   gmshModelGeoMeshSetTransfiniteSurface(302, "Left", cornerTags3, 4, &ierr);
//   gmshModelGeoMeshSetRecombine(2,302, 90, &ierr);

//   int cornerTags[] = {108, 109, 1009, 1008};
//   gmshModelGeoMeshSetTransfiniteSurface(300, "Left", cornerTags, 4, &ierr);
//   gmshModelGeoMeshSetRecombine(2,300, 90, &ierr);


//   gmshModelGeoSynchronize(&ierr);
//   gmshModelMeshGenerate(2, &ierr);
  
//   gmshModelMeshSetSmoothing(2, 301, 100, &ierr);

//   gmshModelGeoSynchronize(&ierr);
//   gmshModelMeshClear(0, 0, &ierr);
//   gmshModelMeshGenerate(2, &ierr);
//  gmshWrite("diamond.msh", &ierr);

  int gui = 1;
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-nopopup")) {
      gui = 0;
      break;
    }
  }
  if (gui) gmshFltkRun(&ierr);

  gmshFinalize(&ierr);
  gmshFree(sBlOutTags);
  return 0;
  // free(outDimTags);
}
