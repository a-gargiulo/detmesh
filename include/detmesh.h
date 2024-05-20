#ifndef DETMESH_H
#define DETMESH_H

#include "converter.h"

int run(int argc, char** argv);

void clean_resources(Point* points, Curve* curves, Surface* surfaces, Node* nodes,
                     size_t* n_entity_blocks, Element* elements, size_t* n_entity_blocks_elements);

#endif // DETMESH_H
