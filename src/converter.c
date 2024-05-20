#include "converter.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "diamond.h"
#include "error.h"
#include "mesh.h"

#define CONVERTER_LINE_BUFFER_SIZE 1024

int read_mesh_structure(const char* file_name, int* mesh_structure, size_t n_mesh_structure) {
  FILE* file = fopen(file_name, "r");
  if (file == NULL) {
    log_error("Could not open mesh structure file!", ERROR_COULD_NOT_OPEN_FILE);
    return ERROR_COULD_NOT_OPEN_FILE;
  }

  for (size_t i = 0; i < n_mesh_structure; i = i + 2) {
    fscanf(file, " %d %d ", &mesh_structure[i], &mesh_structure[i + 1]);
  }
  fclose(file);

  return 0;
}

int read_gmsh(const char* file_name, Point** points, size_t* n_points, Curve** curves,
              size_t* n_curves, Surface** surfaces, size_t* n_surfaces, Node** nodes,
              size_t* n_nodes, size_t* n_entity_blocks, Element** elements, size_t* n_entity_blocks_elements) {
  size_t min_node_tag, max_node_tag;
  size_t n_elements, min_element_tag, max_element_tag;

  char line_buffer[CONVERTER_LINE_BUFFER_SIZE];
  char start_category[50];
  char end_category[50];

  // Read gmsh mesh file
  FILE* file = fopen(file_name, "r");
  if (file == NULL) {
    log_error("The gmsh mesh file could not be opened!", ERROR_COULD_NOT_OPEN_FILE);
    return ERROR_COULD_NOT_OPEN_FILE;
  }

  printf("\nGmsh Reader\n##################\n\n");
  while (fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file) != NULL || !feof(file)) {
    // Find start of mesh section
    while (sscanf(line_buffer, " $%s \n", start_category) != 1) {
      fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
    }

    if (strcmp(start_category, "Entities") == 0) {
      printf("Reading Entities\t-->\t");

      // Read first line
      fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
      sscanf(line_buffer, " %zu %zu %zu 0 \n", n_points, n_curves, n_surfaces);

      // Read the rest
      *points = (Point*)malloc(*n_points * sizeof(Point));
      *curves = (Curve*)malloc(*n_curves * sizeof(Curve));
      *surfaces = (Surface*)malloc(*n_surfaces * sizeof(Surface));
      if (*points == NULL || *curves == NULL || *surfaces == NULL) {
        log_error("Failed to allocate memory for surfaces!", ERROR_NULL_POINTER);
        return ERROR_NULL_POINTER;
      }

      // Points
      for (size_t i = 0; i < *n_points; ++i) {
        fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
        sscanf(line_buffer, " %d %lf %lf %lf \n", &(*points)[i].tag, &(*points)[i].x,
               &(*points)[i].y, &(*points)[i].z);
      }

      // Curves
      for (size_t i = 0; i < *n_curves; ++i) {
        fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
        sscanf(line_buffer, " %d %lf %lf %lf %lf %lf %lf 0 2 %d %d \n", &(*curves)[i].tag,
               &(*curves)[i].min_x, &(*curves)[i].min_y, &(*curves)[i].min_z, &(*curves)[i].max_x,
               &(*curves)[i].max_y, &(*curves)[i].max_z, &(*curves)[i].tags_bounding_points[0],
               &(*curves)[i].tags_bounding_points[1]);
      }

      // Surfaces
      for (size_t i = 0; i < *n_surfaces; ++i) {
        fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
        sscanf(line_buffer, " %d %lf %lf %lf %lf %lf %lf 0 4 %d %d %d %d \n", &(*surfaces)[i].tag,
               &(*surfaces)[i].min_x, &(*surfaces)[i].min_y, &(*surfaces)[i].min_z,
               &(*surfaces)[i].max_x, &(*surfaces)[i].max_y, &(*surfaces)[i].max_z,
               &(*surfaces)[i].tags_bounding_curves[0], &(*surfaces)[i].tags_bounding_curves[1],
               &(*surfaces)[i].tags_bounding_curves[2], &(*surfaces)[i].tags_bounding_curves[3]);
      }

      // Find end of mesh section
      while (sscanf(line_buffer, " $%s ", end_category) != 1) {
        fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
      }

      if (strcmp(end_category, "EndEntities") == 0) {
        printf("Done!\n");
      }
    } else if (strcmp(start_category, "Nodes") == 0) {
      printf("Reading Nodes\t\t-->\t");

      // Read first line
      fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
      sscanf(line_buffer, " %zu %zu %zu %zu ", n_entity_blocks, n_nodes, &min_node_tag,
             &max_node_tag);

      // Read the rest
      *nodes = (Node*)malloc(*n_entity_blocks * sizeof(Node));
      if (*nodes == NULL) {
        log_error("Failed to allocate memory for nodes!", ERROR_NULL_POINTER);
        return ERROR_NULL_POINTER;
      }

      for (size_t i = 0; i < *n_entity_blocks; ++i) {
        fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
        sscanf(line_buffer, " %d %d 0 %zu ", &(*nodes)[i].entity_dim, &(*nodes)[i].entity_tag,
               &(*nodes)[i].n_nodes_in_block);

        (*nodes)[i].node_tags = (size_t*)malloc((*nodes)[i].n_nodes_in_block * sizeof(size_t));
        (*nodes)[i].x = (double*)malloc((*nodes)[i].n_nodes_in_block * sizeof(double));
        (*nodes)[i].y = (double*)malloc((*nodes)[i].n_nodes_in_block * sizeof(double));
        (*nodes)[i].z = (double*)malloc((*nodes)[i].n_nodes_in_block * sizeof(double));

        if ((*nodes)[i].node_tags == NULL || (*nodes)[i].x == NULL || (*nodes)[i].y == NULL ||
            (*nodes)[i].z == NULL) {
          log_error("Failed to allocate memory for node tags, x, y, or z!", ERROR_NULL_POINTER);
          return ERROR_NULL_POINTER;
        }

        for (size_t j = 0; j < (*nodes)[i].n_nodes_in_block; ++j) {
          fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
          sscanf(line_buffer, " %zu ", &(*nodes)[i].node_tags[j]);
        }

        for (size_t k = 0; k < (*nodes)[i].n_nodes_in_block; ++k) {
          fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
          sscanf(line_buffer, " %lf %lf %lf ", &(*nodes)[i].x[k], &(*nodes)[i].y[k],
                 &(*nodes)[i].z[k]);
        }
      }

      // Find end of mesh section
      while (sscanf(line_buffer, " $%s ", end_category) != 1) {
        fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
      }

      if (strcmp(end_category, "EndNodes") == 0) {
        printf("Done!\n");
      }
    } else if (strcmp(start_category, "Elements") == 0) {
      printf("Reading Elements\t-->\t");

      // Read first line
      fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
      sscanf(line_buffer, " %zu %zu %zu %zu ", n_entity_blocks_elements, &n_elements,
             &min_element_tag, &max_element_tag);

      // Read the rest
      *elements = (Element*)malloc(*n_entity_blocks_elements * sizeof(Element));
      if (*elements == NULL) {
        log_error("Failed to allocate memory for elements!", ERROR_NULL_POINTER);
        return ERROR_NULL_POINTER;
      }

      for (size_t i = 0; i < *n_entity_blocks_elements; ++i) {
        fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
        sscanf(line_buffer, " %d %d %d %zu ", &(*elements)[i].entity_dim,
               &(*elements)[i].entity_tag, &(*elements)[i].element_type,
               &(*elements)[i].n_elements_in_block);

        (*elements)[i].element_tags =
            (size_t*)malloc((*elements)[i].n_elements_in_block * sizeof(size_t));
        if ((*elements)[i].element_tags == NULL) {
          log_error("Failed to allocate memory for element tags!", ERROR_NULL_POINTER);
          return ERROR_NULL_POINTER;
        }

        if ((*elements)[i].element_type == 15) {
          (*elements)[i].node_tags =
              (size_t*)malloc((*elements)[i].n_elements_in_block * sizeof(size_t));

          if ((*elements)[i].node_tags == NULL) {
            log_error("Failed to allocate memory for node tags!", ERROR_NULL_POINTER);
            return ERROR_NULL_POINTER;
          }

          for (size_t j = 0; j < (*elements)[i].n_elements_in_block; ++j) {
            fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
            sscanf(line_buffer, " %zu %zu ", &(*elements)[i].element_tags[j],
                   &(*elements)[i].node_tags[j]);
          }
        } else if ((*elements)[i].element_type == 1) {
          (*elements)[i].node_tags =
              (size_t*)malloc(2 * (*elements)[i].n_elements_in_block * sizeof(size_t));

          if ((*elements)[i].node_tags == NULL) {
            log_error("Failed to allocate memory for node tags!", ERROR_NULL_POINTER);
            return ERROR_NULL_POINTER;
          }

          for (size_t j = 0; j < (*elements)[i].n_elements_in_block; ++j) {
            fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
            sscanf(line_buffer, " %zu %zu %zu ", &(*elements)[i].element_tags[j],
                   &(*elements)[i].node_tags[j * 2], &(*elements)[i].node_tags[j * 2 + 1]);
          }
        } else if ((*elements)[i].element_type == 3) {
          (*elements)[i].node_tags =
              (size_t*)malloc(4 * (*elements)[i].n_elements_in_block * sizeof(size_t));

          if ((*elements)[i].node_tags == NULL) {
            log_error("Failed to allocate memory for node tags!", ERROR_NULL_POINTER);
            return ERROR_NULL_POINTER;
          }

          for (size_t j = 0; j < (*elements)[i].n_elements_in_block; ++j) {
            fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
            sscanf(line_buffer, " %zu %zu %zu %zu %zu ", &(*elements)[i].element_tags[j],
                   &(*elements)[i].node_tags[j * 4], &(*elements)[i].node_tags[j * 4 + 1],
                   &(*elements)[i].node_tags[j * 4 + 2], &(*elements)[i].node_tags[j * 4 + 3]);
          }
        }
      }
      while (sscanf(line_buffer, " $%s ", end_category) != 1) {
        fgets(line_buffer, CONVERTER_LINE_BUFFER_SIZE, file);
      }

      if (strcmp(end_category, "EndElements") == 0) {
        printf("Done!\n");
      }
    }
  }

  fclose(file);
  return 0;
}

int ySorter(const void* point1, const void* point2) {
  Point2D* pointA = (Point2D*)point1;
  Point2D* pointB = (Point2D*)point2;
  return (pointA->y > pointB->y) - (pointA->y < pointB->y);
}

void transpose(FluentNode* arr, int numRows, int numCols) {
  FluentNode** transposed = (int**)malloc(numCols * sizeof(FluentNode*));
  for (int i = 0; i < numCols; i++) {
    transposed[i] = (FluentNode*)malloc(numRows * sizeof(FluentNode));
  }

  // Transpose the array
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      transposed[j][i] = arr[j * numRows + i];  // Swap rows and columns
    }
  }

  // Copy the transposed array back to the original array
  for (int i = 0; i < numCols; ++i) {
    for (int j = 0; j < numRows; ++j) {
      arr[j * numCols + i] = transposed[i][j];
    }
  }

  for (int i = 0; i < numCols; i++) {
    free(transposed[i]);
  }
  free(transposed);
}

int is_boundary(const Node* nodes, int block) {
  int boundaries[] = {0, 13, 1, 54, 0, 58, 1, 1007, 0, 135, 1, 133, 0, 107};
  int n_boundaries = 14;

  for (int i = 0; i < n_boundaries; i = i + 2) {
    if (nodes[block].entity_dim == boundaries[i] && nodes[block].entity_tag == boundaries[i + 1])
      return 1;
  }
  return 0;
}

int is_reversed(const Node* nodes, int block) {
  int revBlocks[] = {1,   108, 2,   110, 1,   109, 2,   114, 1,   113, 2,   118, 1,   117, 2,
                     122, 1,   121, 2,   126, 1,   125, 2,   130, 1,   129, 2,   134, 1,   133};
  int n_revBlocks = 30;

  for (int i = 0; i < n_revBlocks; i = i + 2) {
    if (nodes[block].entity_dim == revBlocks[i] && nodes[block].entity_tag == revBlocks[i + 1])
      return 1;
  }
  return 0;
}

int write_fluent(const char* output_file, const Node* nodes, const size_t n_nodes, const size_t n_entity_blocks, const Diamond* diamond, const MeshConfig* mesh_config) {
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

  size_t n_fluent_nodes = n_nodes - 1;  // discard arc center
  FluentNode* f_nodes = (FluentNode*)malloc(n_fluent_nodes * sizeof(FluentNode));
  if (f_nodes == NULL) {
    log_error("Failed to allocate memory for Fluent nodes!", ERROR_NULL_POINTER);
    return ERROR_NULL_POINTER;
  }

  printf("\nFluent Mesh Writer\n##################\n\n");

  int* mesh_structure = (int*)calloc(2 * (n_entity_blocks - 1), sizeof(int));
  if (mesh_structure == NULL) {
    log_error("Failed to allocate memory for mesh structure!", ERROR_NULL_POINTER);
    free(f_nodes);
    return ERROR_NULL_POINTER;
  }

  status = read_mesh_structure("../data/structure.txt", mesh_structure, 2 * (n_entity_blocks - 1));
  if (status != 0) {
    free(mesh_structure);
    free(f_nodes);
    return status;
  }

  // Find mesh dimensions 
  dim_y = 0;
  for (size_t i = 0; i < n_entity_blocks; ++i) {
    if (nodes[i].entity_dim == 2)
      continue;
    else {
      for (size_t j = 0; j < nodes[i].n_nodes_in_block; ++j) {
        if (nodes[i].x[j] == 0.0) {
          dim_y++;
        }
      }
    }
  }
  if (dim_y != 0) {
    dim_x = (n_nodes - 1) / dim_y;
  } else {
    log_error("The y-dimension resulted to zero!", ERROR_DIVISION_BY_ZERO);
    free(f_nodes);
    free(mesh_structure);
    return ERROR_DIVISION_BY_ZERO;
  }
  printf("Mesh Dimensions:\n");
  printf("----------------\n");
  printf("x: %d\n", dim_x);
  printf("y: %d\n", dim_y);
  printf("----------------\n\n");

  // Traverse mesh structure
  for (size_t i = 0; i < 2 * (n_entity_blocks - 1); i = i + 2) {
    // Search block of nodes
    block = 0;
    while (nodes[block].entity_dim != mesh_structure[i] ||
           nodes[block].entity_tag != mesh_structure[i + 1]) {
      block++;
    }

    // Reorder reversed blocks
    Point2D* shuffler = (Point2D*)malloc(nodes[block].n_nodes_in_block * sizeof(Point2D));
    if (shuffler == NULL) {
      log_error("Failed to allocate memory for shuffler!", ERROR_NULL_POINTER);
      free(f_nodes);
      free(mesh_structure);
      return ERROR_NULL_POINTER;
    }
    for (size_t m = 0; m < nodes[block].n_nodes_in_block; ++m) {
      shuffler[m].x = nodes[block].x[m];
      shuffler[m].y = nodes[block].y[m];
      shuffler[m].tag = nodes[block].node_tags[m];
    }

    if (is_reversed(nodes, block)) {
      if (nodes[block].entity_dim == 1) {
        qsort(shuffler, nodes[block].n_nodes_in_block, sizeof(Point2D), ySorter);
        if (reversed_blocks_dim_y == 0) {
          reversed_blocks_dim_y = nodes[block].n_nodes_in_block;
        }
      } else if (nodes[block].entity_dim == 2) {
        for (size_t m = 0; m < nodes[block].n_nodes_in_block / reversed_blocks_dim_y; ++m) {
          qsort(shuffler + m * reversed_blocks_dim_y, reversed_blocks_dim_y, sizeof(Point2D), ySorter);
        }
      }
      for (size_t m = 0; m < nodes[block].n_nodes_in_block; ++m) {
        nodes[block].x[m] = shuffler[m].x;
        nodes[block].y[m] = shuffler[m].y;
        nodes[block].node_tags[m] = shuffler[m].tag;
      }
    }

    // Process each block of nodes sequentially
    if (is_horizontal > 0) {
      for (size_t k = 0; k < nodes[block].n_nodes_in_block; ++k) {
        f_nodes[col_idx * dim_x + row_idx].x = nodes[block].x[k];
        f_nodes[col_idx * dim_x + row_idx].y = nodes[block].y[k];
        f_nodes[col_idx * dim_x + row_idx].tag = col_idx * dim_x + row_idx;
        row_idx++;
      }
    }
    // if NOT horizontal
    else {
      if (nodes[block].entity_dim == 1) {
        col_idx_tmp = col_idx;
        for (size_t k = 0; k < nodes[block].n_nodes_in_block; ++k) {
          f_nodes[col_idx_tmp * dim_x + row_idx].x = nodes[block].x[k];
          f_nodes[col_idx_tmp * dim_x + row_idx].y = nodes[block].y[k];
          f_nodes[col_idx_tmp * dim_x + row_idx].tag = col_idx_tmp * dim_x + row_idx;
          col_idx_tmp++;
        }
        row_idx++;
      } else if (nodes[block].entity_dim == 2) {
        col_idx_tmp = col_idx;
        for (size_t k = 0; k < nodes[block].n_nodes_in_block; ++k) {
          f_nodes[col_idx_tmp * dim_x + row_idx].x = nodes[block].x[k];
          f_nodes[col_idx_tmp * dim_x + row_idx].y = nodes[block].y[k];
          f_nodes[col_idx_tmp * dim_x + row_idx].tag = col_idx_tmp * dim_x + row_idx;

          if (k + 1 != nodes[block].n_nodes_in_block) {
            if (nodes[block].y[k + 1] > nodes[block].y[k] &&
                fabs(nodes[block].x[k + 1] - nodes[block].x[k]) < 1e-4) {
              col_idx_tmp++;
            } else {
              col_idx_tmp = col_idx;
              row_idx++;
            }
          } else {
            col_idx_tmp = col_idx;
            row_idx++;
          }
        }
      }
    }

    if (is_boundary(nodes, block)) {
      row_idx = 0;
      if (is_horizontal > 0)
        col_idx++;
      else
        col_idx = col_idx + nodes[block].n_nodes_in_block;
      n++;
      if (n % 2 == 0) {
        is_horizontal = 1;
      } else
        is_horizontal = -1;
    }
    free(shuffler);
  }

  printf("Task 1: Write node data to 'plot_data.txt' for plotting --> ");
  file_tmp = fopen("plot_data.txt", "w");
  for (size_t i = 0; i < n_fluent_nodes; ++i) {
    fprintf(file_tmp, "%.12lf %.12lf\n", f_nodes[i].x, f_nodes[i].y);
  }
  fclose(file_tmp);
  printf("Done!\n\n");

  int n_fluent_cells = (dim_x - 1) * (dim_y - 1);
  int n_fluent_faces = dim_x * (dim_y - 1) + dim_y * (dim_x - 1);
  int n_fluent_faces_top = dim_x - 1;
  int n_fluent_faces_io = dim_y - 1;

  // find nr of diamond and free surface faces on symmetry boundary
  // boundary points
  int bounding_points[2] = {3, 11};
  int n_bounding_points = 2;
  int bounding_blocks[2] = {0, 0};
  for (int i = 0; i < n_bounding_points; ++i) {
    while (nodes[bounding_blocks[i]].entity_dim != 0 || nodes[bounding_blocks[i]].entity_tag != bounding_points[i]) {
      bounding_blocks[i]++;
    }
  }

  int idx_up = 0;
  int idx_down = 0;
  for (int i = 0; i < dim_x; ++i) {
    if (f_nodes[i].x <= *nodes[bounding_blocks[0]].x && f_nodes[i].y == 0.0) {
      idx_up++;
    } else if (f_nodes[i].x >= *nodes[bounding_blocks[1]].x && f_nodes[i].y == 0.0) {
      idx_down++;
    }
  }

  int n_fluent_faces_up = idx_up - 1;
  int n_fluent_faces_down = idx_down - 1;
  int n_fluent_faces_diamond = dim_x - idx_up - idx_down + 1;

  int n_fluent_faces_interior = n_fluent_faces - 2 * n_fluent_faces_io - n_fluent_faces_top -
                               n_fluent_faces_up - n_fluent_faces_diamond - n_fluent_faces_down;

  int row;
  int col;

  int n1;
  int n2;
  int c1;
  int c2;

  int count;
  file = fopen(output_file, "w");

  if (file == NULL) {
    log_error("Could not open the Fluent mesh file", ERROR_COULD_NOT_OPEN_FILE);
    free(f_nodes);
    free(mesh_structure);
    return ERROR_COULD_NOT_OPEN_FILE;
  }

  fprintf(file, "(0 \"Diamond Mesh:\")\n");
  fprintf(file, "\n");
  fprintf(file, "(0 \"Dimensions:\")\n");
  fprintf(file, "(2 2)\n");
  fprintf(file, "\n");
  fprintf(file, "(12 (0 1 %x 0))\n", n_fluent_cells);
  fprintf(file, "(13 (0 1 %x 0))\n", n_fluent_faces);
  fprintf(file, "(10 (0 1 %x 0 2))\n", n_fluent_nodes);
  fprintf(file, "\n");
  fprintf(file, "(12 (2 1 %x 1 3))\n", n_fluent_cells);
  fprintf(file, "\n");

  // interior -- CHECKED!
  count = 0;
  fprintf(file, "(13 (3 1 %x %x 2)(\n", n_fluent_faces_interior, 2);
  for (size_t i = 0; i < n_fluent_nodes; ++i) {
    row = floor(i / dim_x);
    col = i % dim_x;

    // vertical
    if (col > 0 && col < dim_x - 1 && row < dim_y - 1) {
      count++;
      n1 = i + 1;
      n2 = i + dim_x + 1;
      c1 = (floor((i + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + dim_x) % dim_x;
      c2 = c1 + 1;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }

    // horizontal
    if (row > 0 && row < dim_y - 1 && col < dim_x - 1) {
      count++;
      n1 = i + 2;
      n2 = i + 1;
      c2 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
      c1 = c2 - dim_x + 1;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  printf("Interior faces effective: %d\n", count);
  printf("Interior faces: %d\n", n_fluent_faces_interior);
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // inlet -- CHECKED!
  count = 0;
  fprintf(file, "(13 (4 %x %x 4 2)(\n", n_fluent_faces_interior + 1,
          n_fluent_faces_interior + n_fluent_faces_io);
  for (size_t i = 0; i < n_fluent_nodes; ++i) {
    row = floor(i / dim_x);
    col = i % dim_x;

    if (col == 0 && row < dim_y - 1) {
      count++;
      n1 = i + dim_x + 1;
      n2 = i + 1;
      c1 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  printf("Inlet faces effective: %d\n", count);
  printf("Inlet faces: %d\n", n_fluent_faces_io);
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // outlet -- CHECKED!
  fprintf(file, "(13 (5 %x %x 5 2)(\n", n_fluent_faces_interior + n_fluent_faces_io + 1,
          n_fluent_faces_interior + 2 * n_fluent_faces_io);
  count = 0;
  for (size_t i = 0; i < n_fluent_nodes; ++i) {
    row = floor(i / dim_x);
    col = i % dim_x;

    if (col == dim_x - 1 && row < dim_y - 1) {
      count++;
      n1 = i + 1;
      n2 = i + dim_x + 1;
      c1 = (floor((i + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + dim_x) % dim_x;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  printf("Outlet faces effective: %d\n", count);
  printf("Outlet faces: %d\n", n_fluent_faces_io);
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // symmetry up -- CHECKED!
  count = 0;
  fprintf(file, "(13 (6 %x %x 7 2)(\n", n_fluent_faces_interior + 2 * n_fluent_faces_io + 1,
          n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up);
  for (size_t i = 0; i < n_fluent_nodes; ++i) {
    row = floor(i / dim_x);
    col = i % dim_x;

    if (row == 0 && col < dim_x - 1 && f_nodes[i].x < mesh_config->sUp) {
      count++;
      n1 = i + 1;
      n2 = i + 2;
      c1 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  printf("Symmetry Up faces effective: %d\n", count);
  printf("Symmetry Up faces: %d\n", n_fluent_faces_up);
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // diamond -- CHECKED!
  count = 0;
  fprintf(file, "(13 (7 %x %x 3 2)(\n",
          n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + 1,
          n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond);
  for (size_t i = 0; i < n_fluent_nodes; ++i) {
    row = floor(i / dim_x);
    col = i % dim_x;

    if (row == 0 && col < dim_x - 1 && f_nodes[i].x >= mesh_config->sUp &&
        f_nodes[i].x < mesh_config->sUp + diamond->l) {
      count++;
      n1 = i + 1;
      n2 = i + 2;
      c1 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  printf("Diamond faces effective: %d\n", count);
  printf("Diamond faces: %d\n", n_fluent_faces_diamond);
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // symmetry down -- CHECKED!
  count = 0;
  fprintf(
      file, "(13 (8 %x %x 7 2)(\n",
      n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond + 1,
      n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond +
          n_fluent_faces_down);
  for (size_t i = 0; i < n_fluent_nodes; ++i) {
    row = floor(i / dim_x);
    col = i % dim_x;

    if (row == 0 && col < dim_x - 1 && f_nodes[i].x >= mesh_config->sUp + diamond->l) {
      count++;
      n1 = i + 1;
      n2 = i + 2;
      c1 = (floor((i + 1 + dim_x) / dim_x) - 1) * (dim_x - 1) + (i + 1 + dim_x) % dim_x;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  printf("Symmetry down faces effective: %d\n", count);
  printf("Symmetry down faces: %d\n", n_fluent_faces_down);
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // top wall -- CHECKED!
  count = 0;
  fprintf(file, "(13 (9 %x %x 3 2)(\n",
          n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond +
              n_fluent_faces_down + 1,
          n_fluent_faces_interior + 2 * n_fluent_faces_io + n_fluent_faces_up + n_fluent_faces_diamond +
              n_fluent_faces_down + n_fluent_faces_top);
  for (size_t i = 0; i < n_fluent_nodes; ++i) {
    row = floor(i / dim_x);
    col = i % dim_x;

    if (row == dim_y - 1 && col < dim_x - 1) {
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
  printf("Top faces effective: %d\n", count);
  printf("Top faces: %d\n", n_fluent_faces_top);
  printf("CELLS: %d\n", n_fluent_cells);
  fprintf(file, "(10 (1 1 %x 1 2)(\n", n_fluent_nodes);
  for (size_t i = 0; i < n_fluent_nodes; ++i) {
    fprintf(file, "%.6e %.6e\n", f_nodes[i].x, f_nodes[i].y);
  }
  fprintf(file, "))\n");
  fprintf(file, "\n");

  fclose(file);
  printf("FILE WRITTEN!\n");
  free(f_nodes);
  free(mesh_structure);
  return 0;
}
