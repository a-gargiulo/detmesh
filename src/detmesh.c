#include "detmesh.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "converter.h"
#include "diamond.h"
#include "gmesh.h"
#include "parsing.h"

int run(int argc, char** argv)
{
    // VARIABLES
    int status;

    Diamond diamond;
     
    GMeshConfig gmesh_config;
    GMesh gmesh =
    {
        NULL, 0,
        NULL, 0,
        NULL, 0,
        NULL, 0, 0,
        NULL, 0
    };

    double x_guess[4];

    // PARSER
    status = parse_user_input(argv[argc - 1], &diamond, x_guess, &gmesh_config);
    if (status != 0)
    {
        return status;
    }

    // SOLVE DIAMOND ARC PARAMETERS
    status = calculate_arc_parameters(x_guess, &diamond);
    if (status != 0)
    {
        return status;
    }

    // GENERATE GMSH
    status = mesh_diamond(argc, argv, &diamond, &gmesh_config);
    if (status != 0)
    {
        return status;
    }

    // READ GMSH
    status = read_gmsh("diamond.msh", &gmesh);
    if (status != 0)
    {
        clean_resources(points, curves, surfaces, nodes, &n_entity_blocks, elements,
                        &n_entity_blocks_elements);
        return status;
    }

    // WRITE FLUENT MESH
    status =
        write_fluent("diamond_fluent.msh", nodes, &n_nodes, &n_entity_blocks, &diamond, &gmesh_config);
    if (status != 0)
    {
        clean_resources(points, curves, surfaces, nodes, &n_entity_blocks, elements,
                        &n_entity_blocks_elements);
        return status;
    }

    // CLEANUP
    clean_resources(points, curves, surfaces, nodes, &n_entity_blocks, elements,
                    &n_entity_blocks_elements);

    return 0;
}

void clean_resources(Point* points, Curve* curves, Surface* surfaces, Node* nodes,
                     size_t* n_entity_blocks, Element* elements, size_t* n_entity_blocks_elements)
{

    free(points);

    free(curves);

    free(surfaces);

    if (*n_entity_blocks != 0)
    {
        for (size_t i = 0; i < *n_entity_blocks; ++i)
        {
            free(nodes[i].node_tags);
            free(nodes[i].x);
            free(nodes[i].y);
            free(nodes[i].z);
        }
    }
    free(nodes);

    if (*n_entity_blocks_elements)
    {
        for (size_t i = 0; i < *n_entity_blocks_elements; ++i)
        {
            free(elements[i].element_tags);
            free(elements[i].node_tags);
        }
    }
    free(elements);
}
