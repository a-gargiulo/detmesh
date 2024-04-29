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

  size_t sBlSteps = 70;
  int sBlElementsPerStep[70];
  for (size_t i = 0; i < sBlSteps; ++i) {
    sBlElementsPerStep[i] = 1;
  }

  size_t sBlHeights_n = 70;
  double sBlHeights[70];
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
  gmshFree(sBlOutTags);

  gmshModelGeoSynchronize(&ierr);
  gmshModelMeshGenerate(2, &ierr);
 
  // BL connections, index family 100
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

  // Add top wall (100-) 
  tagCount = 99;
  size_t tConnect_n = 7;  
  for (size_t i = 0; i < tConnect_n; ++i) 
    gmshModelGeoAddLine(i + 100, i + 101, ++tagCount, &ierr);

  // top bl
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

  //bottom curves: c12 p(16 17), c16 p(17 21), c20 p(21 25), c24 p(25 29), c28 p(29, 33), c32 p(33 38), c36 p(38 42), c40 p(42 46), c44 p(46 50), c48 p(50 54), c52 p(54 58)
  // top curves:c107 (110 111), c111 (111 115), c115 (115 119), c119 (119 123), c123 (123 127), c127 (127 131), c131 (131 135)
  // bottom Bl connecto points 16, 17, 25, 29, 42, 46, 54, 58
  // BL connectors index 1000 family
  // gmshModelGeoSynchronize(&ierr);

  int blConnectionPts[] = {16, 110, 17, 111, 25, 115, 29, 119, 42, 123, 46, 127, 54, 131, 58, 135};
  size_t blConnectionPts_n = 16;
  tagCount = 0;
  for (size_t i = 0; i < blConnectionPts_n; i+=2) {
    gmshModelGeoAddLine(blConnectionPts[i], blConnectionPts[i+1], 1000 + tagCount, &ierr);
    tagCount++;
    printf("TAG %d\n", tagCount);
  }

  gmshModelGeoSynchronize(&ierr);
  gmshModelMeshClear(0, 0, &ierr);

  // 7 loops
  int l1Tags[] = {12, 1001, -107, -1000};
  int l2Tags[] = {16, 20, 1002, -111, -1001};
  int l3Tags[] = {24, 1003, -115, -1002};
  int l4Tags[] = {28, 32, 36, 1004, -119, -1003};
  int l5Tags[] = {40, 1005, -123, -1004};
  int l6Tags[] = {44, 48, 1006, -127, -1005};
  int l7Tags[] = {52, 1007, -131, -1006} ;

  size_t l1Tags_n = 4;
  size_t l2Tags_n = 5;
  size_t l3Tags_n = 4;
  size_t l4Tags_n = 6;
  size_t l5Tags_n = 4;
  size_t l6Tags_n = 5;
  size_t l7Tags_n = 4;

  gmshModelGeoAddCurveLoop(l1Tags, l1Tags_n, 1, 0, &ierr); // 0 is for reorient
  gmshModelGeoAddCurveLoop(l2Tags, l2Tags_n, 2, 0, &ierr); // 0 is for reorient
  gmshModelGeoAddCurveLoop(l3Tags, l3Tags_n, 3, 1, &ierr); // 0 is for reorient
  gmshModelGeoAddCurveLoop(l4Tags, l4Tags_n, 4, 0, &ierr); // 0 is for reorient
  gmshModelGeoAddCurveLoop(l5Tags, l5Tags_n, 5, 0, &ierr); // 0 is for reorient
  gmshModelGeoAddCurveLoop(l6Tags, l6Tags_n, 6, 0, &ierr); // 0 is for reorient
  gmshModelGeoAddCurveLoop(l7Tags, l7Tags_n, 7, 0, &ierr); // 0 is for reorient

  int s1_loops[] = {1};
  int s2_loops[] = {2};
  int s3_loops[] = {3};
  int s4_loops[] = {4};
  int s5_loops[] = {5};
  int s6_loops[] = {6};
  int s7_loops[] = {7};

  int s1_corners[] = {16, 17, 111, 110};
  int s2_corners[] = {17, 25, 115, 111};
  int s3_corners[] = {25, 29, 119, 115};
  int s4_corners[] = {29, 42, 123, 119};
  int s5_corners[] = {42, 46, 127, 123};
  int s6_corners[] = {46, 54, 131, 127};
  int s7_corners[] = {54, 58, 135, 131};

  size_t loops_n = 1;

  gmshModelGeoAddPlaneSurface(s1_loops, loops_n, 300, &ierr);
  gmshModelGeoAddPlaneSurface(s2_loops, loops_n, 301, &ierr);
  gmshModelGeoAddPlaneSurface(s3_loops, loops_n, 302, &ierr);
  gmshModelGeoAddPlaneSurface(s4_loops, loops_n, 303, &ierr);
  gmshModelGeoAddPlaneSurface(s5_loops, loops_n, 304, &ierr);
  gmshModelGeoAddPlaneSurface(s6_loops, loops_n, 305, &ierr);
  gmshModelGeoAddPlaneSurface(s7_loops, loops_n, 306, &ierr);

  gmshModelGeoSynchronize(&ierr);

  // Blk1 
  gmshModelGeoMeshSetTransfiniteCurve(1, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(1001, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(100, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(1000, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteSurface(300, "Left", s1_corners, 4, &ierr);
  gmshModelGeoMeshSetRecombine(2,300, 90, &ierr);

  // Blk2 
  gmshModelGeoMeshSetTransfiniteCurve(2, 51, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(3, 50, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(1002, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(101, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteSurface(301, "Left", s2_corners, 4, &ierr);
  gmshModelGeoMeshSetRecombine(2,301, 90, &ierr);
  gmshModelGeoMeshSetSmoothing(2, 301, 1000, &ierr);
  // gmshModelMeshGenerate(2, &ierr);

  // gmshModelGeoMeshSetSmoothing(2, 301, 1000, &ierr);
  

  // Blk3 
  gmshModelGeoMeshSetTransfiniteCurve(4, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(1003, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(102, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteSurface(302, "Left", s3_corners, 4, &ierr);
  gmshModelGeoMeshSetRecombine(2,302, 90, &ierr);
  gmshModelGeoMeshSetSmoothing(2, 302, 1000, &ierr);
  
  // Blk4 
  gmshModelGeoMeshSetTransfiniteCurve(5, 33, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(6, 35, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(7, 34, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(1004, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(103, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteSurface(303, "Left", s4_corners, 4, &ierr);
  gmshModelGeoMeshSetRecombine(2,303, 90, &ierr);
  gmshModelGeoMeshSetSmoothing(2, 303, 1000, &ierr);

  // Blk5 
  gmshModelGeoMeshSetTransfiniteCurve(8, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(1005, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(104, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteSurface(304, "Left", s5_corners, 4, &ierr);
  gmshModelGeoMeshSetRecombine(2,304, 90, &ierr);
  gmshModelGeoMeshSetSmoothing(2, 304, 1000, &ierr);

  // Blk6 
  gmshModelGeoMeshSetTransfiniteCurve(9, 50, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(10, 51, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(1006, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(105, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteSurface(305, "Left", s6_corners, 4, &ierr);
  gmshModelGeoMeshSetRecombine(2,305, 90, &ierr);
  gmshModelGeoMeshSetSmoothing(2, 305, 1000, &ierr);

  // Blk7 
  gmshModelGeoMeshSetTransfiniteCurve(11, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(1007, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteCurve(106, 100, "Progression", 1, &ierr);
  gmshModelGeoMeshSetTransfiniteSurface(306, "Left", s7_corners, 4, &ierr);
  gmshModelGeoMeshSetRecombine(2,306, 90, &ierr);
  gmshModelGeoMeshSetSmoothing(2, 306, 1000, &ierr);

  gmshModelGeoSynchronize(&ierr);
  gmshModelMeshGenerate(2, &ierr);

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
  if (ierr == 0) {
    printf("GREAT SUCCESS\n");
  }
  return 0;
  // free(outDimTags);
}
