#ifndef DETMESH_H
#define DETMESH_H

#include "fmesh.h"
#include "gmesh.h"

int run(int argc, char** argv);

const char* get_opt(const char* option, int argc, char** argv);

void clean_resources(GMesh* gmesh, FMesh* fmesh);

#endif // DETMESH_H
