#include "mesh.h"

#define _USE_MATH_DEFINES

#include <gmshc.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "diamond.h"

void mesh_diamond(int argc, char** argv, Diamond* diamond) {
  int ierr;
  gmshInitialize(argc, argv, 1, 0, &ierr);

  if (gmshIsInitialized(&ierr)) {
    printf("SUCCESSFUL INITIALIZATION OF GMSH\n");
  } else {
    printf("ERROR!\n");
  }

  gmshModelAdd("detmesh", &ierr);
  const double lc = 1e-2;
  gmshModelGeoAddPoint(0, 0, 0, lc, 1, &ierr);
  gmshModelGeoAddPoint(diamond->x1,
                       tan(diamond->alpha * M_PI / 180.0) * diamond->x1, 0, lc,
                       2, &ierr);
  gmshModelGeoAddPoint(diamond->cx, diamond->cy, 0, lc, 3, &ierr);
  gmshModelGeoAddPoint(
      diamond->x2,
      tan(diamond->beta * M_PI / 180.0) * (diamond->l_d - diamond->x2), 0, lc,
      4, &ierr);
  gmshModelGeoAddPoint(diamond->l_d, 0, 0, lc, 5, &ierr);

  gmshModelGeoAddLine(1, 2, 1, &ierr);
  gmshModelGeoAddCircleArc(2, 3, 4, 2, 0, 0, 0, &ierr);
  gmshModelGeoAddLine(4, 5, 3, &ierr);
  gmshWrite("diamond.msh", &ierr);

  gmshModelGeoSynchronize(&ierr);

  int gui = 1;
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-nopopup")) {
      gui = 0;
      break;
    }
  }
  if (gui) gmshFltkRun(&ierr);

  gmshFinalize(&ierr);
}
