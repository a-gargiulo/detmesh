#ifndef MESH_H
#define MESH_H

#include "diamond.h"


typedef struct {
  double sUp, sDown;
  double sBlkUp, sBlkDown;
  double sArcUp, sArcDown;
  double tHeight;
} MeshConfig;

int mesh_diamond(int argc, char** argv, Diamond* diamond, MeshConfig* meshConfig);


#endif // MESH_H
