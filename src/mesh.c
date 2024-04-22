#include "mesh.h"

#define _USE_MATH_DEFINES

#include <gmshc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diamond.h"

void mesh_diamond(int argc, char** argv, Diamond* diamond) {
  double tw_b_up = 0.2; 
  double tw_b_down = 0.2; 
  double t_height = 0.1016; 



  int ierr;
  gmshInitialize(argc, argv, 1, 0, &ierr);

  if (gmshIsInitialized(&ierr)) {
    printf("SUCCESSFUL INITIALIZATION OF GMSH\n");
  } else {
    printf("ERROR!\n");
  }

  gmshModelAdd("detmesh", &ierr);
  const double lc = 5e-3;
  gmshModelGeoAddPoint(0, 0, 0, lc, 1, &ierr);
  gmshModelGeoAddPoint(tw_b_up, 0, 0, lc, 2, &ierr);
  gmshModelGeoAddPoint(tw_b_up + diamond->x1,
                       tan(diamond->alpha * M_PI / 180.0) * diamond->x1, 0, 1e-3,
                       3, &ierr);
  gmshModelGeoAddPoint(tw_b_up + diamond->cx, diamond->cy, 0, lc, 4, &ierr);
  gmshModelGeoAddPoint(
      tw_b_up + diamond->x2,
      tan(diamond->beta * M_PI / 180.0) * (diamond->l_d - diamond->x2), 0, 1e-3,
      5, &ierr);
  gmshModelGeoAddPoint(tw_b_up + diamond->l_d, 0, 0, lc, 6, &ierr);
  gmshModelGeoAddPoint(tw_b_up + diamond->l_d + tw_b_down, 0, 0, lc, 7, &ierr);
  gmshModelGeoAddPoint(0, t_height, 0, lc, 8, &ierr);
  gmshModelGeoAddPoint(tw_b_up + diamond->l_d + tw_b_down, t_height, 0, lc, 9, &ierr);

  gmshModelGeoAddLine(1, 2, 1, &ierr);
  gmshModelGeoAddLine(2, 3, 2, &ierr);
  gmshModelGeoAddCircleArc(3, 4, 5, 3, 0, 0, 0, &ierr);
  gmshModelGeoAddLine(5, 6, 4, &ierr);
  gmshModelGeoAddLine(6, 7, 5, &ierr);
  gmshModelGeoAddLine(8, 9, 6, &ierr);

  gmshModelGeoSynchronize(&ierr);

  double ratio = 1.05;
  double d_w = -1e-6;

  int dimTags[] = {1, 1, 1, 2, 1, 3, 1, 4, 1, 5};  
  size_t dimTags_n = 10;

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
    }
    else {
      heights[i] = heights[i-1] - (-d_w) * pow(ratio, i); 
    }
  }

  int recombine = 1;
  int second = 0;
  int viewIndex=-1;
  

  // gmshModelGeoExtrude(dimTags, dimTags_n, 0, 0.05, 0, &outDimTags, &outDimTags_n, numElements, numElements_n, heights, heights_n, recombine, &ierr);
  gmshModelGeoExtrudeBoundaryLayer(dimTags, dimTags_n, &outDimTags, &outDimTags_n, numElements, numElements_n, heights, heights_n, recombine, second, viewIndex, &ierr);  


  double ratio2 = 1.1;
  double d_w2 = 1e-6;

  int dimTags2[] = {1, 6};  
  size_t dimTags_n2 = 2;

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
    }
    else {
      heights2[i] = heights2[i-1] + d_w2 * pow(ratio2, i); 
    }
  }

  int recombine2= 1;
  int second2= 0;
  int viewIndex2=-1;
  

  gmshModelGeoExtrudeBoundaryLayer(dimTags2, dimTags_n2, &outDimTags2, &outDimTags_n2, numElements2, numElements_n2, heights2, heights_n2, recombine2, second2, viewIndex2, &ierr);





  int* outDimTags_t = NULL;
  size_t outDimTags_n_t;
  // gmshModelGeoCopy(dimTags, dimTags_n, &outDimTags_t, &outDimTags_n_t, &ierr);
  for (size_t i = 0; i < outDimTags_n; ++i) {
    printf("%d\n", outDimTags[i]);
  }

  // printf("HEIGHT = %.8f\n", heights[heights_n]);
  // gmshModelGeoTranslate(outDimTags_t, outDimTags_n_t, 0.0, -heights[heights_n-1], 0.0, &ierr);


  // Extruded curves bottom 7, 11, 15, 19, 23
  // gmshModelGeoAddLine(12, 13, 100, &ierr); // curve 7
  // gmshModelGeoAddLine(13, 17, 101, &ierr); // curve 11
  // gmshModelGeoAddCircleArc(17, 4, 22, 102, 0, 0, 0, &ierr); // curve 15
  // gmshModelGeoAddLine(22, 26, 103, &ierr); // curve 19
  // gmshModelGeoAddLine(26, 30, 104, &ierr); // curve 23 
  // Extruded curves top 27
  // gmshModelGeoAddLine(33, 34, 105, &ierr); // curve 27
  

  gmshModelGeoAddLine(12, 33, 100, &ierr); // curve 7
  gmshModelGeoAddLine(30, 34, 101, &ierr); // curve 7

  int loopTags[] = {7, 11, 15, 19, 23, 101, -27, -100};
  // int loopTags[] = {-101, -23, -19, -15, -11, -7, 100, 27};
  gmshModelGeoAddCurveLoop(loopTags, 8, 1, 0 ,&ierr);
  int wireTags[] = {1};
  gmshModelGeoAddPlaneSurface(wireTags, 1, 1, &ierr);

  gmshModelGeoSynchronize(&ierr);
  gmshModelMeshGenerate(2,&ierr);

  
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
  // free(outDimTags);
}
