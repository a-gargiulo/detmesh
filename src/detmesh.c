#include "detmesh.h"

#include <stdio.h>
#include <stdlib.h>

#include "converter.h"
#include "diamond.h"
#include "mesh.h"
#include "parsing.h"

int run(int argc, char **argv)
{
    int status;

    double x_guess[4];

    Diamond diamond;
    MeshConfig mesh_config;
    
    int n_points, n_curves, n_surfaces;
    int n_entity_blocks, n_entity_blocks_element;
    int numNodes;

    Point *points = NULL;
    Curve *curves = NULL;
    Surface *surfaces = NULL;
    Node *nodes = NULL;
    Element *elements = NULL;

    // PARSE INPUT FILE
    status = parse_user_input(argv[argc - 1], &diamond, x_guess, &mesh_config);
    if (status != 0)
    {
        return status;
    }

    // CALCULATE DIAMOND ARC PARAMETERS
    status = calculate_arc_parameters(x_guess, &diamond);
    if (status != 0)
    {
        return status;
    }

    // GENERATE GMSH
    status = mesh_diamond(argc, argv, &diamond, &mesh_config);
    if (status != 0)
    {
        return status;
    }

    // READ GMSH FILE
    read_gmsh("diamond.msh", &points, &n_points, &curves, &n_curves, &surfaces, &n_surfaces, &nodes, &n_entity_blocks,
             &elements, &n_entity_blocks_element, &numNodes);

    // WRITE FLUENT MESH
    writeFluent("diamond_fluent.msh", nodes, n_entity_blocks, numNodes, &diamond, &mesh_config);

    // CLEANUP
    free(points);
    free(curves);
    free(surfaces);
    for (int i = 0; i < n_entity_blocks; ++i)
    {
        free(nodes[i].nodeTags);
        free(nodes[i].x);
        free(nodes[i].y);
        free(nodes[i].z);
    }
    free(nodes);
    for (int i = 0; i < n_entity_blocks_element; ++i)
    {
        free(elements[i].elementTags);
        free(elements[i].nodeTags);
    }
    free(elements);
    return 0;
}
