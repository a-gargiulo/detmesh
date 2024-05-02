#ifndef CONVERTER_H
#define CONVERTER_H
#include "diamond.h"
#include "mesh.h"
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
  int numNodesInBlock;
  int* nodeTags;
  double* x;
  double* y;
  double* z;
} Node;

typedef struct {
  int entityDim;
  int entityTag;
  int elementType;
  int numElementsInBlock;
  int* elementTags;
  int* nodeTags;
} Element;

typedef struct {
  double x, y;
  int tag;
} FluentNode;


int readGmsh(const char* fileName, Point** points, int* numPoints, Curve** curves, int* numCurves, Surface** surfaces, int* numSurfaces, Node** nodes, int* numEntityBlocks, Element** elements, int* numEntityBlocksElem, int* numNodes);


void print_title();

int writeFluent(const char* outputFile, const Node* nodes, const int numEntityBlocks, const int numNodes, const Diamond* diamond, const MeshConfig* meshConfig);

int xSorter(const void* node1, const void* node2);
int ySorter(const void* node1, const void* node2);

void transpose(FluentNode* arr, int numCols, int numRows);

#endif // CONVERTER_H
