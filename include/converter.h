#ifndef CONVERTER_H
#define CONVERTER_H

#include <stddef.h>


typedef struct {
  int tag;
  double x, y, z;
} Point; 

typedef struct {
  int tag;
  double minX, minY, minZ;
  double maxX, maxY, maxZ;
  int tagsBoundingPoints[2]; // lines
} Curve;

typedef struct {
  int tag;
  double minX, minY, minZ;
  double maxX, maxY, maxZ;
  int tagsBoundingCurves[4]; // quads
} Surface;


typedef struct {
  int entityDim;
  int entityTag;
  size_t numNodesInBlock;
  size_t* nodeTags;
  double* x, *y, *z;
} Node;


void readGmsh(const char* fileName, Point** points, size_t* numPoints, Curve** curves, size_t* numCurves, Surface** surfaces, size_t* numSurfaces, Node** nodes);



#endif // CONVERTER_H
