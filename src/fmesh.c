#include "fmesh.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diamond.h"
#include "error.h"
#include "gmesh.h"

#define CONVERTER_LINE_BUFFER_SIZE 1024
#define CONVERTER_CATEGORY_BUFFER_SIZE 50
#define CONVERTER_N_STRUCTURE_CATEGORIES 3

struct FluentNode
{
    int tag;
    double x, y;
};

struct Point2D
{
    int tag;
    double x, y;
};

struct StructureFileElement
{
    const char* category;
    size_t n_memory;
    int* memory;
};

static StructureFileElement structure_file_elements[CONVERTER_N_STRUCTURE_CATEGORIES] = 
{
    {"OrderEntities", 0, NULL},
    {"OutletEntities", 0, NULL},
    {"ReversedEntities", 0, NULL}
};

// MESH STRUCTURE FILE
int read_structure_file_block(FILE* file, char* line_buffer, size_t n_line_buffer,
                              StructureFileElement* structure_file_element)
{
    fgets(line_buffer, n_line_buffer, file);
    sscanf(line_buffer, " %zu \n", &structure_file_element->n_memory);

    structure_file_element->memory = (int*)malloc(structure_file_element->n_memory * sizeof(int));
    if (structure_file_element->memory == NULL)
    {
        log_error("Could not allocate memory for a structure file element!", ERROR_NULL_POINTER);
        return ERROR_NULL_POINTER;
    }

    for (size_t i = 0; i < structure_file_element->n_memory; i = i + 2)
    {
        fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
        sscanf(line_buffer, " %d %d \n", &structure_file_element->memory[i], &structure_file_element->memory[i + 1]);
    }

    return 0;
}

int read_mesh_structure_file(const char* file_name)
{
    int status;

    char line_buffer[CONVERTER_LINE_BUFFER_SIZE];
    char category_buffer[CONVERTER_CATEGORY_BUFFER_SIZE];

    FILE* file = fopen(file_name, "r");
    if (file == NULL)
    {
        log_error("Mesh structure file could not be opened!", ERROR_COULD_NOT_OPEN_FILE);
        return ERROR_COULD_NOT_OPEN_FILE;
    }

    while (fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file) != NULL || !feof(file))
    {
        for (size_t i = 0; i < CONVERTER_N_STRUCTURE_CATEGORIES; ++i)
        {
            while (sscanf(line_buffer, " $%s \n", category_buffer) != 1)
            {
                fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
            }

            if (strcmp(structure_file_elements[i].category, category_buffer) == 0)
            {
                status = read_structure_file_block(file, line_buffer, CONVERTER_LINE_BUFFER_SIZE,
                                                   &structure_file_elements[i]);
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
    for (size_t i = 0; i < CONVERTER_N_STRUCTURE_CATEGORIES; ++i)
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
    for (size_t i = 0; i < CONVERTER_N_STRUCTURE_CATEGORIES; ++i)
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
        if (nodes[block].entity_dim == outlet_entities[i] && nodes[block].entity_tag == outlet_entities[i + 1])
            return 1;
    }
    return 0;
}

int is_reversed(const Node* nodes, int block, const int* reversed_entities, size_t n_reversed_entities)
{
    for (size_t i = 0; i < n_reversed_entities; i = i + 2)
    {
        if (nodes[block].entity_dim == reversed_entities[i] && nodes[block].entity_tag == reversed_entities[i + 1])
            return 1;
    }
    return 0;
}

int write_fluent(const char* output_file, const Node* nodes, const size_t* n_nodes, const size_t* n_entity_blocks,
                 const Diamond* diamond, const MeshConfig* mesh_config)
{
    int status;

    int dim_x, dim_y;

    int block;
    int n = 0;

    int row_idx = 0;
    int col_idx = 0;
    int col_idx_tmp = 0;
    int is_horizontal = 1;
    int reversed_blocks_dim_y = 0;

    FILE* file;
    FILE* file_tmp;

    size_t n_fluent_nodes = *n_nodes - 1; // discard arc center
    // size_t n_order_entities = 2 * (*n_entity_blocks - 1);

    printf("\nFLUENT MESH WRITER\n##################\n\n");

    FluentNode* f_nodes = (FluentNode*)malloc(n_fluent_nodes * sizeof(FluentNode));
    if (f_nodes == NULL)
    {
        log_error("Failed to allocate memory for Fluent nodes!", ERROR_NULL_POINTER);
        return ERROR_NULL_POINTER;
    }

    status = read_mesh_structure_file("../data/structure.txt");
    if (status != 0)
    {
        for (size_t i = 0; i < CONVERTER_N_STRUCTURE_CATEGORIES; ++i)
        {
            free(structure_file_elements[i].memory);
        }
        free(f_nodes);
        return status;
    }

    // Find mesh dimensions
    dim_y = 0;
    for (size_t i = 0; i < *n_entity_blocks; ++i)
    {
        if (nodes[i].entity_dim == 2)
            continue;
        else
        {
            for (size_t j = 0; j < nodes[i].n_nodes_in_block; ++j)
            {
                if (nodes[i].x[j] == 0.0)
                {
                    dim_y++;
                }
            }
        }
    }
    if (dim_y != 0)
    {
        dim_x = (*n_nodes - 1) / dim_y;
    }
    else
    {
        log_error("The y-dimension resulted to zero!", ERROR_DIVISION_BY_ZERO);
        free(f_nodes);
        for (size_t i = 0; i < CONVERTER_N_STRUCTURE_CATEGORIES; ++i)
        {
            free(structure_file_elements[i].memory);
        }
        return ERROR_DIVISION_BY_ZERO;
    }
    printf("Mesh Dimensions:\n");
    printf("----------------\n");
    printf("x: %d\n", dim_x);
    printf("y: %d\n", dim_y);
    printf("----------------\n\n");

    // Traverse mesh structure
    for (size_t i = 0; i < 2 * (*n_entity_blocks - 1); i = i + 2)
    {
        // Search block of nodes
        block = 0;
        while (nodes[block].entity_dim != get_mesh_structure_element("OrderEntities")[i] ||
               nodes[block].entity_tag != get_mesh_structure_element("OrderEntities")[i + 1])
        {
            block++;
        }

        if (is_reversed(nodes, block, get_mesh_structure_element("ReversedEntities"),
                        get_mesh_structure_element_size("ReversedEntities")))
        {
            // Reorder reversed blocks
            Point2D* shuffler = (Point2D*)malloc(nodes[block].n_nodes_in_block * sizeof(Point2D));
            if (shuffler == NULL)
            {
                log_error("Failed to allocate memory for shuffler!", ERROR_NULL_POINTER);
                free(f_nodes);
                for (size_t i = 0; i < CONVERTER_N_STRUCTURE_CATEGORIES; ++i)
                {
                    free(structure_file_elements[i].memory);
                }
                return ERROR_NULL_POINTER;
            }
            for (size_t m = 0; m < nodes[block].n_nodes_in_block; ++m)
            {
                shuffler[m].x = nodes[block].x[m];
                shuffler[m].y = nodes[block].y[m];
                shuffler[m].tag = nodes[block].node_tags[m];
            }

            if (nodes[block].entity_dim == 1)
            {
                qsort(shuffler, nodes[block].n_nodes_in_block, sizeof(Point2D), y_sorter);
                if (reversed_blocks_dim_y == 0)
                {
                    reversed_blocks_dim_y = nodes[block].n_nodes_in_block;
                }
            }
            else if (nodes[block].entity_dim == 2)
            {
                for (size_t m = 0; m < nodes[block].n_nodes_in_block / reversed_blocks_dim_y; ++m)
                {
                    qsort(shuffler + m * reversed_blocks_dim_y, reversed_blocks_dim_y, sizeof(Point2D), y_sorter);
                }
            }
            for (size_t m = 0; m < nodes[block].n_nodes_in_block; ++m)
            {
                nodes[block].x[m] = shuffler[m].x;
                nodes[block].y[m] = shuffler[m].y;
                nodes[block].node_tags[m] = shuffler[m].tag;
            }
            free(shuffler);
        }

        // Process each block of nodes sequentially
        if (is_horizontal > 0)
        {
            for (size_t k = 0; k < nodes[block].n_nodes_in_block; ++k)
            {
                f_nodes[col_idx * dim_x + row_idx].x = nodes[block].x[k];
                f_nodes[col_idx * dim_x + row_idx].y = nodes[block].y[k];
                f_nodes[col_idx * dim_x + row_idx].tag = col_idx * dim_x + row_idx;
                row_idx++;
            }
        }
        // if NOT horizontal
        else
        {
            if (nodes[block].entity_dim == 1)
            {
                col_idx_tmp = col_idx;
                for (size_t k = 0; k < nodes[block].n_nodes_in_block; ++k)
                {
                    f_nodes[col_idx_tmp * dim_x + row_idx].x = nodes[block].x[k];
                    f_nodes[col_idx_tmp * dim_x + row_idx].y = nodes[block].y[k];
                    f_nodes[col_idx_tmp * dim_x + row_idx].tag = col_idx_tmp * dim_x + row_idx;
                    col_idx_tmp++;
                }
                row_idx++;
            }
            else if (nodes[block].entity_dim == 2)
            {
                col_idx_tmp = col_idx;
                for (size_t k = 0; k < nodes[block].n_nodes_in_block; ++k)
                {
                    f_nodes[col_idx_tmp * dim_x + row_idx].x = nodes[block].x[k];
                    f_nodes[col_idx_tmp * dim_x + row_idx].y = nodes[block].y[k];
                    f_nodes[col_idx_tmp * dim_x + row_idx].tag = col_idx_tmp * dim_x + row_idx;

                    if (k + 1 != nodes[block].n_nodes_in_block)
                    {
                        if (nodes[block].y[k + 1] > nodes[block].y[k] &&
                            fabs(nodes[block].x[k + 1] - nodes[block].x[k]) < 1e-4)
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

        if (is_boundary(nodes, block, get_mesh_structure_element("OutletEntities"),
                        get_mesh_structure_element_size("OutletEntities")))
        {
            row_idx = 0;
            if (is_horizontal > 0)
                col_idx++;
            else
                col_idx = col_idx + nodes[block].n_nodes_in_block;
            n++;
            if (n % 2 == 0)
            {
                is_horizontal = 1;
            }
            else
                is_horizontal = -1;
        }
    }

    printf("Task 1: Write node data to 'plot_data.txt' for plotting --> ");
    file_tmp = fopen("plot_data.txt", "w");
    if (file_tmp == NULL)
    {
        log_error("Coult not open plot_data.txt", ERROR_COULD_NOT_OPEN_FILE);
        free(f_nodes);
        for (size_t i = 0; i < CONVERTER_N_STRUCTURE_CATEGORIES; ++i)
        {
            free(structure_file_elements[i].memory);
        }
        return ERROR_COULD_NOT_OPEN_FILE;
    }
    for (size_t i = 0; i < n_fluent_nodes; ++i)
    {
        fprintf(file_tmp, "%.12lf %.12lf\n", f_nodes[i].x, f_nodes[i].y);
    }
    fclose(file_tmp);
    printf("Done!\n\n");

    // WRITE THE FLUENT MESH FILE
    int n_fluent_cells = (dim_x - 1) * (dim_y - 1);
    int n_fluent_faces = dim_x * (dim_y - 1) + dim_y * (dim_x - 1);
    int n_fluent_faces_top = dim_x - 1;
    int n_fluent_faces_io = dim_y - 1;

    // find nr of diamond and free surface faces on symmetry boundary
    // boundary points
    int bounding_points[2] = {3, 11};
    int n_bounding_points = 2;
    int bounding_blocks[2] = {0, 0};
    for (int i = 0; i < n_bounding_points; ++i)
    {
        while (nodes[bounding_blocks[i]].entity_dim != 0 || nodes[bounding_blocks[i]].entity_tag != bounding_points[i])
        {
            bounding_blocks[i]++;
        }
    }

    int idx_up = 0;
    int idx_down = 0;
    for (int i = 0; i < dim_x; ++i)
    {
        if (f_nodes[i].x <= *nodes[bounding_blocks[0]].x && f_nodes[i].y == 0.0)
        {
            idx_up++;
        }
        else if (f_nodes[i].x >= *nodes[bounding_blocks[1]].x && f_nodes[i].y == 0.0)
        {
            idx_down++;
        }
    }

    int n_fluent_faces_up = idx_up - 1;
    int n_fluent_faces_down = idx_down - 1;
    int n_fluent_faces_diamond = dim_x - idx_up - idx_down + 1;

    int n_fluent_faces_interior = n_fluent_faces - 2 * n_fluent_faces_io - n_fluent_faces_top - n_fluent_faces_up -
                                  n_fluent_faces_diamond - n_fluent_faces_down;

    int row;
    int col;

    int n1;
    int n2;
    int c1;
    int c2;

    int count;
    printf("Task 2: Write %s --> ", output_file);
    file = fopen(output_file, "w");

    if (file == NULL)
    {
        log_error("Could not open the Fluent mesh file", ERROR_COULD_NOT_OPEN_FILE);
        free(f_nodes);
        for (size_t i = 0; i < CONVERTER_N_STRUCTURE_CATEGORIES; ++i)
        {
            free(structure_file_elements[i].memory);
        }
        return ERROR_COULD_NOT_OPEN_FILE;
    }

    fprintf(file, "(0 \"Diamond Mesh:\")\n");
    fprintf(file, "\n");
    fprintf(file, "(0 \"Dimensions:\")\n");
    fprintf(file, "(2 2)\n");
    fprintf(file, "\n");
    fprintf(file, "(12 (0 1 %x 0))\n", n_fluent_cells);
    fprintf(file, "(13 (0 1 %x 0))\n", n_fluent_faces);
    fprintf(file, "(10 (0 1 %x 0 2))\n", (int)n_fluent_nodes);
    fprintf(file, "\n");
    fprintf(file, "(12 (2 1 %x 1 3))\n", n_fluent_cells);
    fprintf(file, "\n");

    // interior -- CHECKED!
    count = 0;
    fprintf(file, "(13 (3 1 %x %x 2)(\n", n_fluent_faces_interior, 2);
    for (size_t i = 0; i < n_fluent_nodes; ++i)
    {
        row = floor(i / dim_x);
        col = i % dim_x;

        // vertical
        if (col > 0 && col < dim_x - 1 && row < dim_y - 1)
        {
            count++;
            n1 = i + 1;
            n2 = i + dim_x + 1;
            c1 = (floor((i + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + dim_x) % dim_x;
            c2 = c1 + 1;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }

        // horizontal
        if (row > 0 && row < dim_y - 1 && col < dim_x - 1)
        {
            count++;
            n1 = i + 2;
            n2 = i + 1;
            c2 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
            c1 = c2 - dim_x + 1;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // inlet -- CHECKED!
    count = 0;
    fprintf(file, "(13 (4 %x %x 4 2)(\n", n_fluent_faces_interior + 1, n_fluent_faces_interior + n_fluent_faces_io);
    for (size_t i = 0; i < n_fluent_nodes; ++i)
    {
        row = floor(i / dim_x);
        col = i % dim_x;

        if (col == 0 && row < dim_y - 1)
        {
            count++;
            n1 = i + dim_x + 1;
            n2 = i + 1;
            c1 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // outlet -- CHECKED!
    fprintf(file, "(13 (5 %x %x 5 2)(\n", n_fluent_faces_interior + n_fluent_faces_io + 1,
            n_fluent_faces_interior + 2 * n_fluent_faces_io);
    count = 0;
    for (size_t i = 0; i < n_fluent_nodes; ++i)
    {
        row = floor(i / dim_x);
        col = i % dim_x;

        if (col == dim_x - 1 && row < dim_y - 1)
        {
            count++;
            n1 = i + 1;
            n2 = i + dim_x + 1;
            c1 = (floor((i + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + dim_x) % dim_x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // symmetry up -- CHECKED!
    count = 0;
    fprintf(file, "(13 (6 %x %x 7 2)(\n", n_fluent_faces_interior + 2 * n_fluent_faces_io + 1,
            n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up);
    for (size_t i = 0; i < n_fluent_nodes; ++i)
    {
        row = floor(i / dim_x);
        col = i % dim_x;

        if (row == 0 && col < dim_x - 1 && f_nodes[i].x < mesh_config->sUp)
        {
            count++;
            n1 = i + 1;
            n2 = i + 2;
            c1 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // diamond -- CHECKED!
    count = 0;
    fprintf(file, "(13 (7 %x %x 3 2)(\n", n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + 1,
            n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond);
    for (size_t i = 0; i < n_fluent_nodes; ++i)
    {
        row = floor(i / dim_x);
        col = i % dim_x;

        if (row == 0 && col < dim_x - 1 && f_nodes[i].x >= mesh_config->sUp &&
            f_nodes[i].x < mesh_config->sUp + diamond->l)
        {
            count++;
            n1 = i + 1;
            n2 = i + 2;
            c1 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // symmetry down -- CHECKED!
    count = 0;
    fprintf(file, "(13 (8 %x %x 7 2)(\n",
            n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond + 1,
            n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond +
                n_fluent_faces_down);
    for (size_t i = 0; i < n_fluent_nodes; ++i)
    {
        row = floor(i / dim_x);
        col = i % dim_x;

        if (row == 0 && col < dim_x - 1 && f_nodes[i].x >= mesh_config->sUp + diamond->l)
        {
            count++;
            n1 = i + 1;
            n2 = i + 2;
            c1 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    // top wall -- CHECKED!
    count = 0;
    fprintf(file, "(13 (9 %x %x 3 2)(\n",
            n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond +
                n_fluent_faces_down + 1,
            n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond +
                n_fluent_faces_down + n_fluent_faces_top);
    for (size_t i = 0; i < n_fluent_nodes; ++i)
    {
        row = floor(i / dim_x);
        col = i % dim_x;

        if (row == dim_y - 1 && col < dim_x - 1)
        {
            count++;
            n1 = i + 2;
            n2 = i + 1;
            c1 = (floor((i + 1) / dim_x) - 1) * (dim_x - 1) + (i + 1) % dim_x;
            c2 = 0;
            fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
        }
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");
    fprintf(file, "(10 (1 1 %x 1 2)(\n", (int)n_fluent_nodes);
    for (size_t i = 0; i < n_fluent_nodes; ++i)
    {
        fprintf(file, "%.6e %.6e\n", f_nodes[i].x, f_nodes[i].y);
    }
    fprintf(file, "))\n");
    fprintf(file, "\n");

    fclose(file);
    printf("Done!\n\n");
    printf("Mesh Stats:\n-----------\n\n");
    printf("Faces: %d\n", n_fluent_faces);
    printf("\t- Interior: %d\n", n_fluent_faces_interior);
    printf("\t- Inlet: %d\n", n_fluent_faces_io);
    printf("\t- Outlet: %d\n", n_fluent_faces_io);
    printf("\t- Symmetry Up: %d\n", n_fluent_faces_up);
    printf("\t- Diamond: %d\n", n_fluent_faces_diamond);
    printf("\t- Symmetry Down: %d\n", n_fluent_faces_down);
    printf("\t- Top Faces: %d\n\n", n_fluent_faces_top);
    printf("Cells: %d\n", n_fluent_cells);
    printf("-----------\n\n");
    free(f_nodes);
    for (size_t i = 0; i < CONVERTER_N_STRUCTURE_CATEGORIES; ++i)
    {
        free(structure_file_elements[i].memory);
    }
    return 0;
}
