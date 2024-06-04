#ifndef FMESH_H
#define FMESH_H

#include <stddef.h>
#include <stdio.h>

#include "diamond.h"
#include "gmesh.h"

typedef struct 
{
    int tag;
    double x, y;
} FluentNode;

typedef struct
{
    int tag;
    double x, y;
} Point2D;

typedef struct
{
    const char* category;
    size_t n_memory;
    int* memory;
} StructureFileElement;

typedef struct
{
  int x, y;
} Dim;

typedef struct
{
  Dim dim;
  size_t  n_nodes;
  size_t n_cells;
  size_t n_faces;
  size_t n_faces_top;
  size_t n_faces_up, n_faces_diamond, n_faces_down;
  size_t n_faces_io;
  size_t n_faces_interior;
  FluentNode* nodes;
} FMesh;

void print_fmesh_stats(const FMesh* fmesh);
int write_fluent_mesh_file(const char* output_file, const Diamond* diamond, const GMeshConfig* gmesh_config, const FMesh* fmesh);

int write_fluent(const char* output_file, const char* structure_file, const Diamond* diamond, const GMeshConfig* gmesh_config, const GMesh* gmesh, FMesh* fmesh);

int y_sorter(const void* node_1, const void* node_2);

int read_mesh_structure_file(const char* file_name);

int read_structure_file_block(FILE* file, char* line_buffer, StructureFileElement* structure_file_element);

int* get_mesh_structure_element(const char* key);

size_t get_mesh_structure_element_size(const char* key);

int is_boundary(const Node* nodes, int block, const int* outlet_entities, size_t n_outlet_entities);

int is_reversed(const Node* nodes, int block, const int* reversed_entities, size_t n_reversed_entities);

void clean_interal_resources(StructureFileElement* structure_file_elements);

int calculate_fmesh_dimensions(const GMesh* gmesh, Dim* dim);

int build_fmesh(const GMesh* gmesh, FMesh* fmesh);

#endif // FMESH_H
