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
  double min_x, max_x;
  double min_y, max_y;
  double min_z, max_z;
  int tags_bounding_points[2];
} Curve;

typedef struct {
  int tag;
  double min_x, max_x;
  double min_y, max_y;
  double min_z, max_z;
  int tags_bounding_curves[4];
} Surface;

typedef struct {
  int entity_dim;
  int entity_tag;
  size_t n_nodes_in_block;
  size_t* node_tags;
  double* x;
  double* y;
  double* z;
} Node;

typedef struct {
  int entity_dim;
  int entity_tag;
  int element_type;
  size_t n_elements_in_block;
  size_t* element_tags;
  size_t* node_tags;
} Element;

typedef struct {
  int tag;
  double x, y;
} FluentNode;

typedef struct {
  int tag;
  double x, y;
} Point2D;

int read_gmsh(const char* file_name, Point** points, size_t* n_points,
              Curve** curves, size_t* n_curves, Surface** surfaces,
              size_t* n_surfaces, Node** nodes, size_t* n_nodes, size_t* n_entity_blocks,
              Element** elements, size_t* n_entity_blocks_elements);

int write_fluent(const char* output_file, const Node* nodes, const size_t n_nodes,
                const size_t n_entity_blocks,
                const Diamond* diamond, const MeshConfig* mesh_config);

int ySorter(const void* node1, const void* node2);

void transpose(FluentNode* arr, int numCols, int numRows);
int read_mesh_structure(const char* file_name, int* mesh_structure,
                       size_t n_mesh_structure);

int is_boundary(const Node* nodes, int block);
int is_reversed(const Node* nodes, int block);
#endif  // CONVERTER_H
