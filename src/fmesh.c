#include "fmesh.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diamond.h"
#include "error.h"
#include "gmesh.h"

#define FMESH_LINE_BUFFER_SIZE 1024
#define FMESH_CATEGORY_BUFFER_SIZE 50
#define FMESH_N_STRUCTURE_CATEGORIES 4
#define FMESH_N_DIAMOND_BOUNDING_ENTITIES 2

static StructureFileElement structure_file_elements[FMESH_N_STRUCTURE_CATEGORIES] = 
{
    {"OrderEntities", 0, NULL},
    {"OutletEntities", 0, NULL},
    {"ReversedEntities", 0, NULL},
    {"DiamondBoundingEntities", 0, NULL}
};

int read_structure_file_block(FILE* file, char* line_buffer,
                              StructureFileElement* structure_file_element)
{
    fgets(line_buffer, FMESH_LINE_BUFFER_SIZE, file);
    sscanf(line_buffer, " %zu \n", &structure_file_element->n_memory);

    structure_file_element->memory = (int*)malloc(structure_file_element->n_memory * sizeof(int));
    if (structure_file_element->memory == NULL)
    {
        log_error("Could not allocate memory for a structure file element!", ERROR_NULL_POINTER);
        return ERROR_NULL_POINTER;
    }

    for (size_t i = 0; i < structure_file_element->n_memory; i = i + 2)
    {
        fgets(line_buffer, FMESH_LINE_BUFFER_SIZE, file);
        sscanf(line_buffer, " %d %d \n", &structure_file_element->memory[i],
               &structure_file_element->memory[i + 1]);
    }

    return 0;
}

int read_mesh_structure_file(const char* file_name)
{
    int status;

    char line_buffer[FMESH_LINE_BUFFER_SIZE];
    char category_buffer[FMESH_CATEGORY_BUFFER_SIZE];

    FILE* file = fopen(file_name, "r");
    if (file == NULL)
    {
        log_error("Mesh structure file could not be opened!", ERROR_COULD_NOT_OPEN_FILE);
        return ERROR_COULD_NOT_OPEN_FILE;
    }

    while (fgets(line_buffer, FMESH_LINE_BUFFER_SIZE, file) != NULL || !feof(file))
    {
        for (size_t i = 0; i < FMESH_N_STRUCTURE_CATEGORIES; ++i)
        {
            while (sscanf(line_buffer, " $%s \n", category_buffer) != 1)
            {
                fgets(line_buffer, FMESH_LINE_BUFFER_SIZE, file);
            }

            if (strcmp(structure_file_elements[i].category, category_buffer) == 0)
            {
                status = read_structure_file_block(file, line_buffer, &structure_file_elements[i]);
                if (status != 0)
                {
                    return status;
                }
            }
        }
    }

    fclose(file);

    return 0;
}

int* get_mesh_structure_element(const char* key)
{
    for (size_t i = 0; i < FMESH_N_STRUCTURE_CATEGORIES; ++i)
    {
        if (strcmp(structure_file_elements[i].category, key) == 0)
        {
            return structure_file_elements[i].memory;
        }
    }

    return NULL;
}

size_t get_mesh_structure_element_size(const char* key)
{
    for (size_t i = 0; i < FMESH_N_STRUCTURE_CATEGORIES; ++i)
    {
        if (strcmp(structure_file_elements[i].category, key) == 0)
        {
            return structure_file_elements[i].n_memory;
        }
    }

    return 0;
}

int y_sorter(const void* node_1, const void* node_2)
{
    Point2D* point_A = (Point2D*)node_1;
    Point2D* point_B = (Point2D*)node_2;
    return (point_A->y > point_B->y) - (point_A->y < point_B->y);
}

int is_boundary(const Node* nodes, int block, const int* outlet_entities, size_t n_outlet_entities)
{
    for (size_t i = 0; i < n_outlet_entities; i = i + 2)
    {
        if (nodes[block].entity_dim == outlet_entities[i] &&
            nodes[block].entity_tag == outlet_entities[i + 1])
            return 1;
    }
    return 0;
}

int is_reversed(const Node* nodes, int block, const int* reversed_entities,
                size_t n_reversed_entities)
{
    for (size_t i = 0; i < n_reversed_entities; i = i + 2)
    {
        if (nodes[block].entity_dim == reversed_entities[i] &&
            nodes[block].entity_tag == reversed_entities[i + 1])
            return 1;
    }
    return 0;
}

void clean_internal_resources(StructureFileElement* structure_file_elements)
{
    for (size_t i = 0; i < FMESH_N_STRUCTURE_CATEGORIES; ++i)
    {
        free(structure_file_elements[i].memory);
    }
}

int calculate_fmesh_dimensions(const GMesh* gmesh, Dim* dim)
{
    dim->y = 0;
    for (size_t i = 0; i < gmesh->n_entity_blocks; ++i)
    {
        if (gmesh->nodes[i].entity_dim == 2)
            continue;
        else
        {
            for (size_t j = 0; j < gmesh->nodes[i].n_nodes_in_block; ++j)
            {
                if (gmesh->nodes[i].x[j] == 0.0)
                {
                    dim->y++;
                }
            }
        }
    }

    if (dim->y == 0)
    {
        log_error("The y-dimension resulted to zero!", ERROR_DIVISION_BY_ZERO);
        return ERROR_DIVISION_BY_ZERO;
    }
    
    dim->x = (gmesh->n_nodes - 1) / dim->y;
    return 0;
}

int build_fmesh(const GMesh* gmesh, FMesh* fmesh)
{
    int block;
    int n = 0;
    int is_horizontal = 1;
    int reversed_blocks_dim_y = 0;
    int row_idx = 0;
    int col_idx = 0;
    int col_idx_tmp = 0;
    int diamond_bounding_entities[FMESH_N_DIAMOND_BOUNDING_ENTITIES] = {0};
    int idx_up = 0;
    int idx_down = 0;
    Point2D* shuffler;

    // Traverse mesh structure
    for (size_t i = 0; i < 2 * (gmesh->n_entity_blocks - 1); i = i + 2)
    {
        // Find node block 
        block = 0;
        while (gmesh->nodes[block].entity_dim != get_mesh_structure_element("OrderEntities")[i] ||
               gmesh->nodes[block].entity_tag != get_mesh_structure_element("OrderEntities")[i + 1])
        {
            block++;
        }

        if (is_reversed(gmesh->nodes, block, get_mesh_structure_element("ReversedEntities"),
                        get_mesh_structure_element_size("ReversedEntities")))
        {
            // Reorder reversed blocks
            shuffler = (Point2D*)malloc(gmesh->nodes[block].n_nodes_in_block * sizeof(Point2D));
            if (shuffler == NULL)
            {
                log_error("Failed to allocate memory for shuffler!", ERROR_NULL_POINTER);
                return ERROR_NULL_POINTER;
            }

            for (size_t m = 0; m < gmesh->nodes[block].n_nodes_in_block; ++m)
            {
                shuffler[m].x = gmesh->nodes[block].x[m];
                shuffler[m].y = gmesh->nodes[block].y[m];
                shuffler[m].tag = gmesh->nodes[block].node_tags[m];
            }

            if (gmesh->nodes[block].entity_dim == 1)
            {
                qsort(shuffler, gmesh->nodes[block].n_nodes_in_block, sizeof(Point2D), y_sorter);
                if (reversed_blocks_dim_y == 0)
                {
                    reversed_blocks_dim_y = gmesh->nodes[block].n_nodes_in_block;
                }
            }
            else if (gmesh->nodes[block].entity_dim == 2)
            {
                for (size_t m = 0; m < gmesh->nodes[block].n_nodes_in_block / reversed_blocks_dim_y;
                     ++m)
                {
                    qsort(shuffler + m * reversed_blocks_dim_y, reversed_blocks_dim_y,
                          sizeof(Point2D), y_sorter);
                }
            }
            for (size_t m = 0; m < gmesh->nodes[block].n_nodes_in_block; ++m)
            {
                gmesh->nodes[block].x[m] = shuffler[m].x;
                gmesh->nodes[block].y[m] = shuffler[m].y;
                gmesh->nodes[block].node_tags[m] = shuffler[m].tag;
            }
            free(shuffler);
        }

        // Process each block of nodes sequentially
        if (is_horizontal > 0)
        {
            for (size_t k = 0; k < gmesh->nodes[block].n_nodes_in_block; ++k)
            {
                fmesh->nodes[col_idx * fmesh->dim.x + row_idx].x = gmesh->nodes[block].x[k];
                fmesh->nodes[col_idx * fmesh->dim.x + row_idx].y = gmesh->nodes[block].y[k];
                fmesh->nodes[col_idx * fmesh->dim.x + row_idx].tag = col_idx * fmesh->dim.x + row_idx;
                row_idx++;
            }
        }
        else
        {
            if (gmesh->nodes[block].entity_dim == 1)
            {
                col_idx_tmp = col_idx;
                for (size_t k = 0; k < gmesh->nodes[block].n_nodes_in_block; ++k)
                {
                    fmesh->nodes[col_idx_tmp * fmesh->dim.x + row_idx].x = gmesh->nodes[block].x[k];
                    fmesh->nodes[col_idx_tmp * fmesh->dim.x + row_idx].y = gmesh->nodes[block].y[k];
                    fmesh->nodes[col_idx_tmp * fmesh->dim.x + row_idx].tag = col_idx_tmp * fmesh->dim.x + row_idx;
                    col_idx_tmp++;
                }
                row_idx++;
            }
            else if (gmesh->nodes[block].entity_dim == 2)
            {
                col_idx_tmp = col_idx;
                for (size_t k = 0; k < gmesh->nodes[block].n_nodes_in_block; ++k)
                {
                    fmesh->nodes[col_idx_tmp * fmesh->dim.x + row_idx].x = gmesh->nodes[block].x[k];
                    fmesh->nodes[col_idx_tmp * fmesh->dim.x + row_idx].y = gmesh->nodes[block].y[k];
                    fmesh->nodes[col_idx_tmp * fmesh->dim.x + row_idx].tag = col_idx_tmp * fmesh->dim.x + row_idx;

                    if (k + 1 != gmesh->nodes[block].n_nodes_in_block)
                    {
                        if (gmesh->nodes[block].y[k + 1] > gmesh->nodes[block].y[k] &&
                            fabs(gmesh->nodes[block].x[k + 1] - gmesh->nodes[block].x[k]) < 1e-4)
                        {
                            col_idx_tmp++;
                        }
                        else
                        {
                            col_idx_tmp = col_idx;
                            row_idx++;
                        }
                    }
                    else
                    {
                        col_idx_tmp = col_idx;
                        row_idx++;
                    }
                }
            }
        }

        if (is_boundary(gmesh->nodes, block, get_mesh_structure_element("OutletEntities"),
                        get_mesh_structure_element_size("OutletEntities")))
        {
            row_idx = 0;
            if (is_horizontal > 0)
                col_idx++;
            else
                col_idx = col_idx + gmesh->nodes[block].n_nodes_in_block;
            n++;
            if (n % 2 == 0)
            {
                is_horizontal = 1;
            }
            else
                is_horizontal = -1;
        }
    }



    fmesh->n_cells = (fmesh->dim.x - 1) * (fmesh->dim.y - 1);
    fmesh->n_faces = fmesh->dim.x * (fmesh->dim.y - 1) + fmesh->dim.y * (fmesh->dim.x - 1);
    fmesh->n_faces_top = fmesh->dim.x - 1;

    for (size_t i = 0; i < get_mesh_structure_element_size("DiamondBoundingEntities") / 2; ++i)
    {
        while (gmesh->nodes[diamond_bounding_entities[i]].entity_dim != 0 ||
               gmesh->nodes[diamond_bounding_entities[i]].entity_tag != get_mesh_structure_element("DiamondBoundingEntities")[2 * i + 1])
        {
            diamond_bounding_entities[i]++;
        }
    }

    for (int i = 0; i < fmesh->dim.x; ++i)
    {
        if (fmesh->nodes[i].x <= *gmesh->nodes[diamond_bounding_entities[0]].x && fmesh->nodes[i].y == 0.0)
        {
            idx_up++;
        }
        else if (fmesh->nodes[i].x >= *gmesh->nodes[diamond_bounding_entities[1]].x && fmesh->nodes[i].y == 0.0)
        {
            idx_down++;
        }
    }

    fmesh->n_faces_up = idx_up - 1;
    fmesh->n_faces_down = idx_down - 1;
    fmesh->n_faces_diamond = fmesh->dim.x - idx_up - idx_down + 1;
    fmesh->n_faces_io = fmesh->dim.y - 1;
    fmesh->n_faces_interior = fmesh->n_faces - 2 * fmesh->n_faces_io - fmesh->n_faces_top -
                                  fmesh->n_faces_up - fmesh->n_faces_diamond - fmesh->n_faces_down;


    return 0;

}

int write_fluent_mesh_file(const char* output_file, const Diamond* diamond, const GMeshConfig* gmesh_config, const FMesh* fmesh)
{

    int row;
    int col;

    int n1;
    int n2;
    int c1;
    int c2;

    printf("Writing Fluent mesh to %s --> ", output_file);
    FILE* file = fopen(output_file, "w");
    if (file == NULL)
    {
        log_error("Could not open the Fluent mesh file", ERROR_COULD_NOT_OPEN_FILE);
        return ERROR_COULD_NOT_OPEN_FILE;
    }

    fprintf(file, "(0 \"Diamond Mesh:\")\n");
    fprintf(file, "\n");
    fprintf(file, "(0 \"Dimensions:\")\n");
    fprintf(file, "(2 2)\n");
    fprintf(file, "\n");
    fprintf(file, "(12 (0 1 %x 0))\n", (int)fmesh->n_cells);
    fprintf(file, "(13 (0 1 %x 0))\n", (int)fmesh->n_faces);
    fprintf(file, "(10 (0 1 %x 0 2))\n", (int)fmesh->n_nodes);
    fprintf(file, "\n");
    fprintf(file, "(12 (2 1 %x 1 3))\n", (int)fmesh->n_cells);
    fprintf(file, "\n");

    // INTERIOR 
    fprintf(file, "(13 (3 1 %x %x 2)(\n", (int)fmesh->n_faces_interior, 2);
    for (size_t i = 0; i < fmesh->n_nodes; ++i)
    {
        row = floor(i / fmesh->dim.x);
        col = i % fmesh->dim.x;

        // vertical
        if (col > 0 && col < fmesh->dim.x - 1 && row < fmesh->dim.y - 1)
        {
            n1 = i + 1;
            n2 = i + fmesh->dim.x + 1;
            c1 = (floor((i + fmesh->dim.x) / fmesh->dim.x) - 1) * (fmesh->dim.x - 1) + (i + fmesh->dim.x) % fmesh->dim.x;
            c2 = c1 + 1;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }

        // horizontal
        if (row > 0 && row < fmesh->dim.y - 1 && col < fmesh->dim.x - 1)
        {
            n1 = i + 2;
            n2 = i + 1;
            c2 = (floor((i + 1 + fmesh->dim.x) / fmesh->dim.x) - 1) * (fmesh->dim.x - 1) + (i + 1 + fmesh->dim.x) % fmesh->dim.x;
            c1 = c2 - fmesh->dim.x + 1;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // INLET
    fprintf(file, "(13 (4 %x %x 4 2)(\n", (int)fmesh->n_faces_interior + 1,
            (int)(fmesh->n_faces_interior + fmesh->n_faces_io));
    for (size_t i = 0; i < fmesh->n_nodes; ++i)
    {
        row = floor(i / fmesh->dim.x);
        col = i % fmesh->dim.x;

        if (col == 0 && row < fmesh->dim.y - 1)
        {
            n1 = i + fmesh->dim.x + 1;
            n2 = i + 1;
            c1 = (floor((i + 1 + fmesh->dim.x) / fmesh->dim.x) - 1) * (fmesh->dim.x - 1) + (i + 1 + fmesh->dim.x) % fmesh->dim.x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // OUTLET
    fprintf(file, "(13 (5 %x %x 5 2)(\n", (int)(fmesh->n_faces_interior + fmesh->n_faces_io + 1),
            (int)(fmesh->n_faces_interior + 2 * fmesh->n_faces_io));
    for (size_t i = 0; i < fmesh->n_nodes; ++i)
    {
        row = floor(i / fmesh->dim.x);
        col = i % fmesh->dim.x;

        if (col == fmesh->dim.x - 1 && row < fmesh->dim.y - 1)
        {
            n1 = i + 1;
            n2 = i + fmesh->dim.x + 1;
            c1 = (floor((i + fmesh->dim.x) / fmesh->dim.x) - 1) * (fmesh->dim.x - 1) + (i + fmesh->dim.x) % fmesh->dim.x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // SYMMETRY UPSTREAM
    fprintf(file, "(13 (6 %x %x 7 2)(\n", (int)(fmesh->n_faces_interior + 2 * fmesh->n_faces_io + 1),
            (int)(fmesh->n_faces_interior + 2 * fmesh->n_faces_io + fmesh->n_faces_up));
    for (size_t i = 0; i < fmesh->n_nodes; ++i)
    {
        row = floor(i / fmesh->dim.x);
        col = i % fmesh->dim.x;

        if (row == 0 && col < fmesh->dim.x - 1 && fmesh->nodes[i].x < gmesh_config->sUp)
        {
            n1 = i + 1;
            n2 = i + 2;
            c1 = (floor((i + 1 + fmesh->dim.x) / fmesh->dim.x) - 1) * (fmesh->dim.x - 1) + (i + 1 + fmesh->dim.x) % fmesh->dim.x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // DIAMOND
    fprintf(file, "(13 (7 %x %x 3 2)(\n",
            (int)(fmesh->n_faces_interior + 2 * fmesh->n_faces_io + fmesh->n_faces_up + 1),
            (int)(fmesh->n_faces_interior + 2 * fmesh->n_faces_io + fmesh->n_faces_up +
                fmesh->n_faces_diamond));
    for (size_t i = 0; i < fmesh->n_nodes; ++i)
    {
        row = floor(i / fmesh->dim.x);
        col = i % fmesh->dim.x;

        if (row == 0 && col < fmesh->dim.x - 1 && fmesh->nodes[i].x >= gmesh_config->sUp &&
            fmesh->nodes[i].x < gmesh_config->sUp + diamond->l)
        {
            n1 = i + 1;
            n2 = i + 2;
            c1 = (floor((i + 1 + fmesh->dim.x) / fmesh->dim.x) - 1) * (fmesh->dim.x - 1) + (i + 1 + fmesh->dim.x) % fmesh->dim.x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // SYMMETRY DOWN
    fprintf(file, "(13 (8 %x %x 7 2)(\n",
            (int)(fmesh->n_faces_interior + 2 * fmesh->n_faces_io + fmesh->n_faces_up +
                fmesh->n_faces_diamond + 1),
            (int)(fmesh->n_faces_interior + 2 * fmesh->n_faces_io + fmesh->n_faces_up +
                fmesh->n_faces_diamond + fmesh->n_faces_down));
    for (size_t i = 0; i < fmesh->n_nodes; ++i)
    {
        row = floor(i / fmesh->dim.x);
        col = i % fmesh->dim.x;

        if (row == 0 && col < fmesh->dim.x - 1 && fmesh->nodes[i].x >= gmesh_config->sUp + diamond->l)
        {
            n1 = i + 1;
            n2 = i + 2;
            c1 = (floor((i + 1 + fmesh->dim.x) / fmesh->dim.x) - 1) * (fmesh->dim.x - 1) + (i + 1 + fmesh->dim.x) % fmesh->dim.x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // TOP WALL
    fprintf(file, "(13 (9 %x %x 3 2)(\n",
            (int)(fmesh->n_faces_interior + 2 * fmesh->n_faces_io + fmesh->n_faces_up +
                fmesh->n_faces_diamond + fmesh->n_faces_down + 1),
            (int)(fmesh->n_faces_interior + 2 * fmesh->n_faces_io + fmesh->n_faces_up +
                fmesh->n_faces_diamond + fmesh->n_faces_down + fmesh->n_faces_top));
    for (size_t i = 0; i < fmesh->n_nodes; ++i)
    {
        row = floor(i / fmesh->dim.x);
        col = i % fmesh->dim.x;

        if (row == fmesh->dim.y - 1 && col < fmesh->dim.x - 1)
        {
            n1 = i + 2;
            n2 = i + 1;
            c1 = (floor((i + 1) / fmesh->dim.x) - 1) * (fmesh->dim.x - 1) + (i + 1) % fmesh->dim.x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");
    
    // NODES
    fprintf(file, "(10 (1 1 %x 1 2)(\n", (int)fmesh->n_nodes);
    for (size_t i = 0; i < fmesh->n_nodes; ++i)
    {
        fprintf(file, "%.6e %.6e\n", fmesh->nodes[i].x, fmesh->nodes[i].y);
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    fclose(file);

    printf("Done!\n\n");

    return 0;
}


void print_fmesh_stats(const FMesh* fmesh)
{

    printf("MESH STATISTICS\n###############\n\n");
    printf("Dimensions:\n");
    printf("----------------\n");
    printf("x: %d\n", fmesh->dim.x);
    printf("y: %d\n", fmesh->dim.y);
    printf("----------------\n\n");
    printf("Structure:\n---------\n\n");
    printf("Nodes: %zu\n", fmesh->n_nodes);
    printf("Cells: %zu\n", fmesh->n_cells);
    printf("Faces: %zu\n", fmesh->n_faces);
    printf("\t- Interior: %zu\n", fmesh->n_faces_interior);
    printf("\t- Inlet: %zu\n", fmesh->n_faces_io);
    printf("\t- Outlet: %zu\n", fmesh->n_faces_io);
    printf("\t- Symmetry Up: %zu\n", fmesh->n_faces_up);
    printf("\t- Diamond: %zu\n", fmesh->n_faces_diamond);
    printf("\t- Symmetry Down: %zu\n", fmesh->n_faces_down);
    printf("\t- Top Faces: %zu\n\n", fmesh->n_faces_top);
    printf("-----------\n\n");
}


int write_fluent(const char* output_file, const Diamond* diamond, const GMeshConfig* gmesh_config,
                 const GMesh* gmesh, FMesh* fmesh)
{
    int status;

    printf("\nFLUENT MESH WRITER\n##################\n\n");

    fmesh->n_nodes = gmesh->n_nodes - 1; // discard arc center
    fmesh->nodes = (FluentNode*)malloc(fmesh->n_nodes * sizeof(FluentNode));
    if (fmesh->nodes == NULL)
    {
        log_error("Failed to allocate memory for Fluent nodes!", ERROR_NULL_POINTER);
        return ERROR_NULL_POINTER;
    }

    status = calculate_fmesh_dimensions(gmesh, &fmesh->dim); 
    if (status != 0)
    {
        clean_internal_resources(structure_file_elements);
        return status;
    }

    status = read_mesh_structure_file("../data/structure.txt");
    if (status != 0)
    {
        clean_internal_resources(structure_file_elements);
        return status;
    }

    status = build_fmesh(gmesh, fmesh); 
    if (status != 0)
    {
        clean_internal_resources(structure_file_elements);
        return status;
    }


    status = write_fluent_mesh_file(output_file, diamond, gmesh_config, fmesh);
    if (status != 0)
    {
        clean_internal_resources(structure_file_elements);
        return status;
    }

    print_fmesh_stats(fmesh);

    clean_internal_resources(structure_file_elements);
    
    return 0;
}
