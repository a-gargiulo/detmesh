#ifndef GMESH_H
#define GMESH_H

#include <stddef.h>
#include <stdio.h>

#include "diamond.h"


typedef struct
{
    int tag;
    double x, y, z;
} Point;

typedef struct
{
    int tag;
    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;
    int tags_bounding_points[2];
} Curve;

typedef struct 
{
    int tag;
    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;
    int tags_bounding_curves[4];
} Surface;

typedef struct
{
    int entity_dim;
    int entity_tag;
    size_t n_nodes_in_block;
    size_t* node_tags;
    double* x;
    double* y;
    double* z;
} Node;

typedef struct
{
    int entity_dim;
    int entity_tag;
    int element_type;
    size_t n_elements_in_block;
    size_t* element_tags;
    size_t* node_tags;
} Element;

typedef struct
{
    Point* points;
    size_t n_points;

    Curve* curves;
    size_t n_curves;

    Surface* surfaces;
    size_t n_surfaces;

    Node* nodes;
    size_t n_nodes, n_entity_blocks;

    Element* elements;
    size_t n_entity_blocks_elements;

} GMesh;

typedef struct
{
    double sUp, sDown;
    double sBlkUp, sBlkDown;
    double sArcUp, sArcDown;
    double tHeight;
} GMeshConfig;

typedef struct
{
    const char* category;
    int (*read_gmsh_type)(FILE*, char*, GMesh*);
} GMeshCategory;


int generate_gmsh(int argc, char** argv, Diamond* diamond, GMeshConfig* gmesh_config);

int read_gmsh(const char* file_name, GMesh* gmesh);

int read_gmsh_entities(FILE* file, char* line_buffer, GMesh* gmesh);

int read_gmsh_nodes(FILE* file, char* line_buffer, GMesh* gmesh);

int read_gmsh_elements(FILE* file, char* line_buffer, GMesh* gmesh);

#endif // GMESH_H
