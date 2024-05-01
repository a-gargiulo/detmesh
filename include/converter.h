#ifndef CONVERTER_H
#define CONVERTER_H

#include <stddef.h>


typedef struct {
  int tag;
  double x; 
  double y; 
  double z;
} Point; 

typedef struct {
  int tag;
  double minX; 
  double minY;
  double minZ;
  double maxX;
  double maxY;
  double maxZ;
  int tagsBoundingPoints[2]; // lines
} Curve;

typedef struct {
  int tag;
  double minX;
  double minY;
  double minZ;
  double maxX;
  double maxY;
  double maxZ;
  int tagsBoundingCurves[4]; // quads
} Surface;


typedef struct {
  int entityDim;
  int entityTag;
  double* x;
  double* y;
  double* z;
  size_t* nodeTags;
  size_t numNodesInBlock;
} Node;

typedef struct {
  int entityDim;
  int entityTag;
  int elementType;
  size_t numElementsInBlock;
  size_t* elementTags;
  size_t* nodeTags;
} Element;

void readGmsh(const char* fileName, Point** points, size_t* numPoints, Curve** curves, size_t* numCurves, Surface** surfaces, size_t* numSurfaces, Node** nodes, size_t* numEntityBlocks, Element** elements, size_t* numEntityBlocksElem);



#endif // CONVERTER_H
