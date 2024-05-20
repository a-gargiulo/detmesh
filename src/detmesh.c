#include "detmesh.h"

#include <stdio.h>
#include <stdlib.h>

#include "converter.h"
#include "diamond.h"
#include "mesh.h"
#include "parsing.h"

int run(int argc, char** argv)
{
    // VARIABLES
    int status;

    Diamond diamond;
    MeshConfig mesh_config;
    double x_guess[4];

    Point* points = NULL;
    Curve* curves = NULL;
    Surface* surfaces = NULL;
    Node* nodes = NULL;
    Element* elements = NULL;
    size_t n_points, n_curves, n_surfaces, n_nodes;
    size_t n_entity_blocks, n_entity_blocks_elements;

    // PARSER
    status = parse_user_input(argv[argc - 1], &diamond, x_guess, &mesh_config);
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
    status = mesh_diamond(argc, argv, &diamond, &mesh_config);
    if (status != 0)
    {
        return status;
    }

    // READ GMSH
    status =
        read_gmsh("diamond.msh", &points, &n_points, &curves, &n_curves, &surfaces, &n_surfaces,
                  &nodes, &n_nodes, &n_entity_blocks, &elements, &n_entity_blocks_elements);
    if (status != 0)
    {
        clean_resources(points, curves, surfaces, nodes, &n_entity_blocks, elements,
                        &n_entity_blocks_elements);
        return status;
    }

    // WRITE FLUENT MESH
    status =
        write_fluent("diamond_fluent.msh", nodes, n_nodes, n_entity_blocks, &diamond, &mesh_config);
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

    for (size_t i = 0; i < *n_entity_blocks; ++i)
    {
        free(nodes[i].node_tags);
        free(nodes[i].x);
        free(nodes[i].y);
        free(nodes[i].z);
    }
    free(nodes);

    for (size_t i = 0; i < *n_entity_blocks_elements; ++i)
    {
        free(elements[i].element_tags);
        free(elements[i].node_tags);
    }
    free(elements);
}
