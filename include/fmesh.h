#ifndef FMESH_H
#define FMESH_H

#include <stddef.h>
#include <stdio.h>

#include "diamond.h"
#include "gmesh.h"

typedef struct FluentNode FluentNode;
typedef struct Point2D Point2D;
typedef struct StructureFileElement StructureFileElement;


int write_fluent(const char* output_file, const Diamond* diamond, const GMeshConfig* gmesh_config, const GMesh* gmesh);

int y_sorter(const void* node_1, const void* node_2);

int read_mesh_structure_file(const char* file_name);

int read_structure_file_block(FILE* file, char* line_buffer, size_t n_line_buffer,
                              StructureFileElement* structure_file_element);

int* get_mesh_structure_element(const char* key);

size_t get_mesh_structure_element_size(const char* key);

int is_boundary(const Node* nodes, int block, const int* outlet_entities, size_t n_outlet_entities);

int is_reversed(const Node* nodes, int block, const int* reversed_entities, size_t n_reversed_entities);

#endif // FMESH_H
