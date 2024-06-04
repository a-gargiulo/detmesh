#include "detmesh.h"

#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include "diamond.h"
#include "fmesh.h"
#include "gmesh.h"
#include "parsing.h"
#include "error.h"

int run(int argc, char** argv)
{
    int status;

    Diamond diamond = {0};
    double x_guess[4] = {0};

    GMeshConfig gmesh_config = {0};
    GMesh gmesh = {0};

    FMesh fmesh = {0};

    
    status = parse_user_input(get_opt("-input", argc, argv), &diamond, x_guess, &gmesh_config);
    if (status != 0)
    {
        clean_resources(&gmesh, &fmesh);
        return status;
    }

    // SOLVE DIAMOND ARC PARAMETERS
    status = calculate_arc_parameters(x_guess, &diamond);
    if (status != 0)
    {
        clean_resources(&gmesh, &fmesh);
        return status;
    }

    // GENERATE GMSH
    status = generate_gmsh(argc, argv, &diamond, &gmesh_config);
    if (status != 0)
    {
        clean_resources(&gmesh, &fmesh);
        return status;
    }

    // READ GMSH
    status = read_gmsh(get_opt("-gmesh", argc, argv), &gmesh);
    if (status != 0)
    {
        clean_resources(&gmesh, &fmesh);
        return status;
    }

    // WRITE FLUENT MESH
    status = write_fluent(get_opt("-fluent", argc, argv), get_opt("-structure", argc, argv), &diamond, &gmesh_config, &gmesh, &fmesh);
    if (status != 0)
    {
        clean_resources(&gmesh, &fmesh);
        return status;
    }

    // CLEANUP
    clean_resources(&gmesh, &fmesh);

    return 0;
}

const char* get_opt(const char* option, int argc, char** argv)
{

    for (int i = 0; i < argc; ++i)
    {
        if (strcmp(argv[i], option) == 0)
            return argv[i+1];
    }

    printf("DETMESH INFO: %s option not found. Using default value.\n", option);
    return NULL;

}


void clean_resources(GMesh* gmesh, FMesh* fmesh)
{

    free(gmesh->points);

    free(gmesh->curves);

    free(gmesh->surfaces);

    if (gmesh->n_entity_blocks != 0)
    {
        for (size_t i = 0; i < gmesh->n_entity_blocks; ++i)
        {
            free(gmesh->nodes[i].node_tags);
            free(gmesh->nodes[i].x);
            free(gmesh->nodes[i].y);
            free(gmesh->nodes[i].z);
        }
    }
    free(gmesh->nodes);

    if (gmesh->n_entity_blocks_elements)
    {
        for (size_t i = 0; i < gmesh->n_entity_blocks_elements; ++i)
        {
            free(gmesh->elements[i].element_tags);
            free(gmesh->elements[i].node_tags);
        }
    }
    free(gmesh->elements);

    free(fmesh->nodes);
}

