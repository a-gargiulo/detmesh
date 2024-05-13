#ifndef CONVERTER_H
#define CONVERTER_H
#include <stddef.h>
#include "diamond.h"
#include "mesh.h"

typedef struct {
  int tag;
  double x, y, z; 
} Point; 

typedef struct {
  int tag;
  double minX, maxX; 
  double minY, maxY; 
  double minZ, maxZ; 
  int tagsBoundingPoints[2]; // lines
} Curve;

typedef struct {
  int tag;
  double minX, maxX; 
  double minY, maxY; 
  double minZ, maxZ; 
  int tagsBoundingCurves[4]; // quads
} Surface;


typedef struct {
  int entityDim;
  int entityTag;
  size_t numNodesInBlock;
  size_t* nodeTags;
  double* x;
  double* y;
  double* z;
} Node;

typedef struct {
  int entityDim;
  int entityTag;
  int elementType;
  size_t numElementsInBlock;
  size_t* elementTags;
  size_t* nodeTags;
} Element;

typedef struct {
  double x, y;
  int tag;
} FluentNode;

typedef struct {
  double x, y;
  int tag;
} Point2D;


int read_gmsh(const char* file_name, Point** points, size_t* n_points, Curve** curves, size_t* n_curves, Surface** surfaces, size_t* n_surfaces, Node** nodes, size_t* n_entity_blocks, Element** elements, size_t* n_entity_blocks_element, size_t* n_nodes);


int writeFluent(const char* outputFile, const Node* nodes, const int numEntityBlocks, const int numNodes, const Diamond* diamond, const MeshConfig* meshConfig);

int ySorter(const void* node1, const void* node2);

void transpose(FluentNode* arr, int numCols, int numRows);
void readMeshStructure(const char* fileName, int* meshStructure,int n_meshStructure);

int isBoundary(const Node* nodes, int block);
int isReversed(const Node* nodes, int block);
#endif // CONVERTER_H
