#include "converter.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "diamond.h"
#include "mesh.h"
#define LINE_BUFFER_SIZE 1024

void readMeshStructure(const char *fileName, int *meshStructure,
                       int n_meshStructure) {
  FILE *file;

  file = fopen(fileName, "r");

  if (file == NULL) {
    fprintf(stderr, "ERROR: Structure file not found");
    return;
  }

  for (int i = 0; i < n_meshStructure; i = i + 2) {
    fscanf(file, " %d %d ", &meshStructure[i], &meshStructure[i + 1]);
  }

  fclose(file);
}

void print_title(void) {
  const char *title =
      "\n"
      "READING GMSH\n"
      "------------\n"
      "\n";
  printf("%s", title);
}

int readGmsh(const char *fileName, Point **points, int *numPoints,
             Curve **curves, int *numCurves, Surface **surfaces,
             int *numSurfaces, Node **nodes, int *numEntityBlocks,
             Element **elements, int *numEntityBlocksElem, int *numNodes) {
  int minNodeTag, maxNodeTag;
  int numElements, minElementTag, maxElementTag;
  FILE *file;
  char line[LINE_BUFFER_SIZE];
  char startCategory[50];
  char endCategory[50];

  print_title();

  printf("Opening file... ");
  file = fopen(fileName, "r");

  if (file == NULL) {
    fprintf(stderr, "ERROR: Gmsh file could not open!\n");
    return 1;
  } else {
    printf("Successful!\n");
  }

  while ((fgets(line, LINE_BUFFER_SIZE, file) != NULL) || (!feof(file))) {
    while (sscanf(line, " $%s \n", startCategory) != 1) {
      fgets(line, LINE_BUFFER_SIZE, file);
    }

    if (strcmp(startCategory, "Entities") == 0) {
      printf("Reading ENTITIES... ");

      // first line
      fgets(line, LINE_BUFFER_SIZE, file);
      sscanf(line, " %d %d %d 0 \n", numPoints, numCurves, numSurfaces);
      *points = (Point *)malloc((*numPoints) * sizeof(Point));
      *curves = (Curve *)malloc((*numCurves) * sizeof(Curve));
      *surfaces = (Surface *)malloc((*numSurfaces) * sizeof(Surface));

      if (*points == NULL || *curves == NULL || *surfaces == NULL) {
        fprintf(stderr,
                "ERROR: Failed to allocate memory for Points, Curves, or "
                "Surfaces!\n");
        return 1;
      }

      // points
      for (int i = 0; i < *numPoints; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %lf %lf %lf \n", &(*points)[i].tag, &(*points)[i].x,
               &(*points)[i].y, &(*points)[i].z);
      }

      // curves
      for (int i = 0; i < *numCurves; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %lf %lf %lf %lf %lf %lf 0 2 %d %d \n",
               &(*curves)[i].tag, &(*curves)[i].minX, &(*curves)[i].minY,
               &(*curves)[i].minZ, &(*curves)[i].maxX, &(*curves)[i].maxY,
               &(*curves)[i].maxZ, &(*curves)[i].tagsBoundingPoints[0],
               &(*curves)[i].tagsBoundingPoints[1]);
      }

      // surfaces
      for (int i = 0; i < *numSurfaces; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %lf %lf %lf %lf %lf %lf 0 4 %d %d %d %d \n",
               &(*surfaces)[i].tag, &(*surfaces)[i].minX, &(*surfaces)[i].minY,
               &(*surfaces)[i].minZ, &(*surfaces)[i].maxX, &(*surfaces)[i].maxY,
               &(*surfaces)[i].maxZ, &(*surfaces)[i].tagsBoundingCurves[0],
               &(*surfaces)[i].tagsBoundingCurves[1],
               &(*surfaces)[i].tagsBoundingCurves[2],
               &(*surfaces)[i].tagsBoundingCurves[3]);
      }

      while (sscanf(line, " $%s ", endCategory) != 1) {
        fgets(line, LINE_BUFFER_SIZE, file);
      }

      if (strcmp(endCategory, "EndEntities") == 0) {
        printf("Reading completed!\n");
      }
    }

    else if (strcmp(startCategory, "Nodes") == 0) {
      printf("Reading NODES... ");

      // first line
      fgets(line, LINE_BUFFER_SIZE, file);
      sscanf(line, " %d %d %d %d ", numEntityBlocks, numNodes, &minNodeTag,
             &maxNodeTag);
      *nodes = (Node *)malloc((*numEntityBlocks) * sizeof(Node));

      if ((*nodes) == NULL) {
        fprintf(stderr, "ERROR: Failed to allocate memory for Nodes!\n");
        return 1;
      }

      // entity blocks
      for (int i = 0; i < *numEntityBlocks; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %d 0 %d ", &(*nodes)[i].entityDim,
               &(*nodes)[i].entityTag, &(*nodes)[i].numNodesInBlock);
        (*nodes)[i].nodeTags =
            (int *)malloc((*nodes)[i].numNodesInBlock * sizeof(int));
        (*nodes)[i].x =
            (double *)malloc((*nodes)[i].numNodesInBlock * sizeof(double));
        (*nodes)[i].y =
            (double *)malloc((*nodes)[i].numNodesInBlock * sizeof(double));
        (*nodes)[i].z =
            (double *)malloc((*nodes)[i].numNodesInBlock * sizeof(double));

        if ((*nodes)[i].nodeTags == NULL || (*nodes)[i].x == NULL ||
            (*nodes)[i].y == NULL || (*nodes)[i].z == NULL) {
          fprintf(stderr,
                  "ERROR: Failed to allocate memory for node tags or x, y, z "
                  "coordinates!\n");
          return 1;
        }

        for (int j = 0; j < (*nodes)[i].numNodesInBlock; ++j) {
          fgets(line, LINE_BUFFER_SIZE, file);
          sscanf(line, " %d ", &(*nodes)[i].nodeTags[j]);
        }

        for (int k = 0; k < (*nodes)[i].numNodesInBlock; ++k) {
          fgets(line, LINE_BUFFER_SIZE, file);
          sscanf(line, " %lf %lf %lf ", &(*nodes)[i].x[k], &(*nodes)[i].y[k],
                 &(*nodes)[i].z[k]);
        }
      }

      while (sscanf(line, " $%s ", endCategory) != 1) {
        fgets(line, LINE_BUFFER_SIZE, file);
      }

      if (strcmp(endCategory, "EndNodes") == 0) {
        printf("Reading completed!\n");
      }
    } else if (strcmp(startCategory, "Elements") == 0) {
      printf("Reading ELEMENTS... ");

      // first line
      fgets(line, LINE_BUFFER_SIZE, file);
      sscanf(line, "%d %d %d %d", numEntityBlocksElem, &numElements,
             &minElementTag, &maxElementTag);
      *elements = (Element *)malloc(*numEntityBlocksElem * sizeof(Element));

      if ((*elements) == NULL) {
        fprintf(stderr, "ERROR: Failed to allocate memory for Elements!\n");
        return 1;
      }

      // Element blocks
      for (int i = 0; i < *numEntityBlocksElem; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %d %d %d ", &(*elements)[i].entityDim,
               &(*elements)[i].entityTag, &(*elements)[i].elementType,
               &(*elements)[i].numElementsInBlock);
        (*elements)[i].elementTags =
            (int *)malloc((*elements)[i].numElementsInBlock * sizeof(int));

        if ((*elements)[i].elementTags == NULL) {
          fprintf(stderr,
                  "ERROR: Failed to allocate memory for element tags!\n");
          return 1;
        }

        if ((*elements)[i].elementType == 15) {
          (*elements)[i].nodeTags =
              (int *)malloc((*elements)[i].numElementsInBlock * sizeof(int));

          if ((*elements)[i].nodeTags == NULL) {
            fprintf(stderr,
                    "ERROR: Failed to allocate memory for node tags!\n");
            return 1;
          }

          for (int j = 0; j < (*elements)[i].numElementsInBlock; ++j) {
            fgets(line, LINE_BUFFER_SIZE, file);
            sscanf(line, " %d %d ", &(*elements)[i].elementTags[j],
                   &(*elements)[i].nodeTags[j]);
          }
        } else if ((*elements)[i].elementType == 1) {
          (*elements)[i].nodeTags = (int *)malloc(
              2 * (*elements)[i].numElementsInBlock * sizeof(int));

          if ((*elements)[i].nodeTags == NULL) {
            fprintf(stderr,
                    "ERROR: Failed to allocate memory for node tags!\n");
            return 1;
          }

          for (int j = 0; j < (*elements)[i].numElementsInBlock; ++j) {
            fgets(line, LINE_BUFFER_SIZE, file);
            sscanf(line, " %d %d %d ", &(*elements)[i].elementTags[j],
                   &(*elements)[i].nodeTags[j * 2],
                   &(*elements)[i].nodeTags[j * 2 + 1]);
          }
        } else if ((*elements)[i].elementType == 3) {
          (*elements)[i].nodeTags = (int *)malloc(
              4 * (*elements)[i].numElementsInBlock * sizeof(int));

          if ((*elements)[i].nodeTags == NULL) {
            fprintf(stderr,
                    "ERROR: Failed to allocate memory for node tags!\n");
            return 1;
          }

          for (int j = 0; j < (*elements)[i].numElementsInBlock; ++j) {
            fgets(line, LINE_BUFFER_SIZE, file);
            sscanf(line, "%d %d %d %d %d ", &(*elements)[i].elementTags[j],
                   &(*elements)[i].nodeTags[j * 4],
                   &(*elements)[i].nodeTags[j * 4 + 1],
                   &(*elements)[i].nodeTags[j * 4 + 2],
                   &(*elements)[i].nodeTags[j * 4 + 3]);
          }
        }
      }
      while (sscanf(line, " $%s ", endCategory) != 1) {
        fgets(line, LINE_BUFFER_SIZE, file);
      }

      if (strcmp(endCategory, "EndElements") == 0) {
        printf("Reading completed!\n");
      }
    }
  }

  printf("Closing file... ");
  int ss = fclose(file);
  if (ss == 0)
    printf("Closed!\n");
  else {
    fprintf(stderr, "ERROR: Could not close file!\n");
    return 1;
  }
  return 0;
}

int ySorter(const void *point1, const void *point2) {
  Point2D *pointA = (Point2D *)point1;
  Point2D *pointB = (Point2D *)point2;
  return (pointA->y > pointB->y) - (pointA->y < pointB->y);
}

void transpose(FluentNode *arr, int numRows, int numCols) {
  FluentNode **transposed = (int **)malloc(numCols * sizeof(FluentNode *));
  for (int i = 0; i < numCols; i++) {
    transposed[i] = (FluentNode *)malloc(numRows * sizeof(FluentNode));
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

int isBoundary(const Node *nodes, int block) {
  int boundaries[] = {0, 13, 1, 54, 0, 58, 1, 1007, 0, 135, 1, 133, 0, 107};
  int n_boundaries = 14;

  for (int i = 0; i < n_boundaries; i = i + 2) {
    if (nodes[block].entityDim == boundaries[i] &&
        nodes[block].entityTag == boundaries[i + 1])
      return 1;
  }
  return 0;
}

int isReversed(const Node *nodes, int block) {
  int revBlocks[] = {1, 108, 2, 110, 1, 109, 2, 114, 1, 113,
                     2, 118, 1, 117, 2, 122, 1, 121, 2, 126,
                     1, 125, 2, 130, 1, 129, 2, 134, 1, 133};
  int n_revBlocks = 30;

  for (int i = 0; i < n_revBlocks; i = i + 2) {
    if (nodes[block].entityDim == revBlocks[i] &&
        nodes[block].entityTag == revBlocks[i + 1])
      return 1;
  }
  return 0;
}

int writeFluent(const char *outputFile, const Node *nodes,
                const int numEntityBlocks, const int numNodes,
                const Diamond *diamond, const MeshConfig *meshConfig) {
  FILE *file;
  int numFluentNodes = numNodes - 1;  // discard arc center
  FluentNode *fNodes =
      (FluentNode *)malloc(numFluentNodes * sizeof(FluentNode));
  if (fNodes == NULL) {
    fprintf(stderr, "Error: Memory allocation failed\n");
    free(fNodes);
    return 1;
  }

  printf("\nWRITING FLUENT\n--------------\n");
  printf("\n# Manipulate Nodes\n##################\n");

  int *meshStructure = (int *)calloc(2 * (numEntityBlocks - 1), sizeof(int));

  readMeshStructure("../data/structure.txt", meshStructure,
                    2 * (numEntityBlocks - 1));

  // i = col_idx * n_rows + row_idx

  int block;
  int n = 0;

  // FIND MESH DIMS
  int dimY = 0;
  for (int i = 0; i < numEntityBlocks; ++i) {
    if (nodes[i].entityDim == 2)
      continue;
    else {
      for (int j = 0; j < nodes[i].numNodesInBlock; ++j) {
        if (nodes[i].x[j] == 0.0) {
          dimY++;
        }
      }
    }
  }
  int dimX;
  if (dimY != 0) {
    dimX = (numNodes - 1) / dimY;
  } else {
    dimX = 0;
    printf("ERROR!");
  }
  printf("DIMENSION Y: %d\n", dimY);
  printf("DIMENSION X: %d\n", dimX);

  int row_idx = 0;
  int col_idx = 0;
  int col_idx_tmp = 0;

  int isHorizontal = 1;
  int revBlockDimY = 0;
  // Traverse mesh structure
  for (int i = 0; i < 2 * (numEntityBlocks - 1); i = i + 2) {
    // Search block of nodes
    block = 0;
    while (nodes[block].entityDim != meshStructure[i] ||
           nodes[block].entityTag != meshStructure[i + 1]) {
      block++;
    }

    // reorder reversed blocks
    Point2D *shuffler =
        (Point2D *)malloc(nodes[block].numNodesInBlock * sizeof(Point2D));
    for (int m = 0; m < nodes[block].numNodesInBlock; ++m) {
      shuffler[m].x = nodes[block].x[m];
      shuffler[m].y = nodes[block].y[m];
      shuffler[m].tag = nodes[block].nodeTags[m];
    }

    if (isReversed(nodes, block)) {
      if (nodes[block].entityDim == 1) {
        qsort(shuffler, nodes[block].numNodesInBlock, sizeof(Point2D), ySorter);
        if (revBlockDimY == 0) {
          revBlockDimY = nodes[block].numNodesInBlock;
        }
      } else if (nodes[block].entityDim == 2) {
        for (int m = 0; m < nodes[block].numNodesInBlock / revBlockDimY; ++m) {
          qsort(shuffler + m * revBlockDimY, revBlockDimY, sizeof(Point2D),
                ySorter);
        }
      }
      for (int m = 0; m < nodes[block].numNodesInBlock; ++m) {
        nodes[block].x[m] = shuffler[m].x;
        nodes[block].y[m] = shuffler[m].y;
        nodes[block].nodeTags[m] = shuffler[m].tag;
      }
    }

    // Process each block of nodes sequentially
    if (isHorizontal > 0) {
      for (int k = 0; k < nodes[block].numNodesInBlock; ++k) {
        fNodes[col_idx * dimX + row_idx].x = nodes[block].x[k];
        fNodes[col_idx * dimX + row_idx].y = nodes[block].y[k];
        fNodes[col_idx * dimX + row_idx].tag = col_idx * dimX + row_idx;
        row_idx++;
      }
    }
    // if NOT horizontal
    else {
      if (nodes[block].entityDim == 1) {
        col_idx_tmp = col_idx;
        for (int k = 0; k < nodes[block].numNodesInBlock; ++k) {
          fNodes[col_idx_tmp * dimX + row_idx].x = nodes[block].x[k];
          fNodes[col_idx_tmp * dimX + row_idx].y = nodes[block].y[k];
          fNodes[col_idx_tmp * dimX + row_idx].tag =
              col_idx_tmp * dimX + row_idx;
          col_idx_tmp++;
        }
        row_idx++;
      } else if (nodes[block].entityDim == 2) {
        col_idx_tmp = col_idx;
        for (int k = 0; k < nodes[block].numNodesInBlock; ++k) {
          fNodes[col_idx_tmp * dimX + row_idx].x = nodes[block].x[k];
          fNodes[col_idx_tmp * dimX + row_idx].y = nodes[block].y[k];
          fNodes[col_idx_tmp * dimX + row_idx].tag =
              col_idx_tmp * dimX + row_idx;

          if (k + 1 != nodes[block].numNodesInBlock) {
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

    if (isBoundary(nodes, block)) {
      row_idx = 0;
      if (isHorizontal > 0)
        col_idx++;
      else
        col_idx = col_idx + nodes[block].numNodesInBlock;
      n++;
      if (n % 2 == 0) {
        isHorizontal = 1;
      } else
        isHorizontal = -1;
    }
    free(shuffler);
  }

  printf("WRITING plot_data.txt... ");
  FILE *tmpFile;
  tmpFile = fopen("plot_data.txt", "w");
  for (int i = 0; i < numFluentNodes; ++i) {
    fprintf(tmpFile, "%.12lf %.12lf\n", fNodes[i].x, fNodes[i].y);
  }
  fclose(tmpFile);
  printf("Done!\n");

  int numFluentCells = (dimX - 1) * (dimY - 1);
  int numFluentFaces = dimX * (dimY - 1) + dimY * (dimX - 1);
  int numFluentFacesTop = dimX - 1;
  int numFluentFacesIO = dimY - 1;

  // find nr of diamond and free surface faces on symmetry boundary
  // boundary points
  int bndPts[2] = {3, 11};
  int n_bndPts = 2;
  int bndBlocks[2] = {0, 0};
  for (int i= 0; i < n_bndPts; ++i) {
    while (nodes[bndBlocks[i]].entityDim != 0 ||
           nodes[bndBlocks[i]].entityTag != bndPts[i]) {
      bndBlocks[i]++;
    }
  }
 
  int idx_up = 0;
  int idx_down = 0;
  for (int i = 0; i < dimX; ++i) {
    if (fNodes[i].x <= *nodes[bndBlocks[0]].x && fNodes[i].y == 0.0) {
      idx_up++;
    } else if (fNodes[i].x >= *nodes[bndBlocks[1]].x && fNodes[i].y == 0.0) {
      idx_down++;
    }
  }

  int numFluentFacesUp = idx_up-1;
  int numFluentFacesDown = idx_down-1;
  int numFluentFacesDiamond = dimX - idx_up - idx_down + 1;

  int numFluentFacesInterior = numFluentFaces - 2 * numFluentFacesIO -
                               numFluentFacesTop - numFluentFacesUp -
                               numFluentFacesDiamond - numFluentFacesDown;

  int row;
  int col;

  int n1;
  int n2;
  int c1;
  int c2;

  file = fopen(outputFile, "w");

  if (file == NULL) {
    fprintf(stderr, "ERROR: Output file could not be opened!\n");
    return 1;
  }

  fprintf(file, "(0 \"Diamond Mesh:\")\n");
  fprintf(file, "\n");
  fprintf(file, "(0 \"Dimensions:\")\n");
  fprintf(file, "(2 2)\n");
  fprintf(file, "\n");
  fprintf(file, "(12 (0 1 %x 0))\n", numFluentCells);
  fprintf(file, "(13 (0 1 %x 0))\n", numFluentFaces);
  fprintf(file, "(10 (0 1 %x 0 2))\n", numFluentNodes);
  fprintf(file, "\n");
  fprintf(file, "(12 (2 1 %x 1 3))\n", numFluentCells);
  fprintf(file, "\n");

  // interior
  fprintf(file, "(13 (3 1 %x %x 2)(\n", numFluentFacesInterior, 2);
  for (int i = 0; i < numFluentNodes; ++i) {
    row = floor(i / dimX);
    col = i % dimX;

    if (col > 0 && col < dimX - 1 && row < dimY - 1) {
      n1 = i + 1;
      n2 = i + dimX + 1;
      c1 = (floor((i + dimX) / dimX) - 1) * (dimX - 1) + (i + dimX) % dimX;
      c2 = c1 + 1;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }

    if (row > 0 && row < dimY - 1 && col < dimX - 1) {
      n1 = i + 2;
      n2 = i + 1;
      c2 = (floor((i + 1 + dimX) / dimX) - 1) * (dimX - 1) +
           (i + 1 + dimX) % dimX;
      c1 = c2 - dimX + 1;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // inlet
  fprintf(file, "(13 (4 %x %x 4 2)(\n", numFluentFacesInterior + 1,
          numFluentFacesInterior + numFluentFacesIO);
  for (int i = 0; i < numFluentNodes; ++i) {
    row = floor(i / dimX);
    col = i % dimX;

    if (col == 0 && row < dimY - 1) {
      n1 = i + dimX + 1;
      n2 = i + 1;
      c1 = (floor((i + 1 + dimX) / dimX) - 1) * (dimX - 1) +
           (i + 1 + dimX) % dimX;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // outlet
  fprintf(file, "(13 (5 %x %x 5 2)(\n",
          numFluentFacesInterior + numFluentFacesIO + 1,
          numFluentFacesInterior + 2 * numFluentFacesIO);
  for (int i = 0; i < numFluentNodes; ++i) {
    row = floor(i / dimX);
    col = i % dimX;

    if (col == dimX - 1 && row < dimY - 1) {
      n1 = i + 1;
      n2 = i + dimX + 1;
      c1 = (floor((i + dimX) / dimX) - 1) * (dimX - 1) + (i + dimX) % dimX;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // symmetry up
  fprintf(file, "(13 (6 %x %x 7 2)(\n",
          numFluentFacesInterior + 2 * numFluentFacesIO + 1,
          numFluentFacesInterior + 2 * numFluentFacesIO + numFluentFacesUp);
  for (int i = 0; i < numFluentNodes; ++i) {
    row = floor(i / dimX);
    col = i % dimX;

    if (row == 0 && col < dimX - 1 && fNodes[i].x < meshConfig->sUp) {
      n1 = i + 1;
      n2 = i + 2;
      c1 = (floor((i + 1 + dimX) / dimX) - 1) * (dimX - 1) +
           (i + 1 + dimX) % dimX;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // diamond
  fprintf(file, "(13 (7 %x %x 3 2)(\n",
          numFluentFacesInterior + 2 * numFluentFacesIO + numFluentFacesUp + 1,
          numFluentFacesInterior + 2 * numFluentFacesIO + numFluentFacesUp +
              numFluentFacesDiamond);
  for (int i = 0; i < numFluentNodes; ++i) {
    row = floor(i / dimX);
    col = i % dimX;

    if (row == 0 && col < dimX - 1 && fNodes[i].x >= meshConfig->sUp &&
        fNodes[i].x < meshConfig->sUp + diamond->l) {
      n1 = i + 1;
      n2 = i + 2;
      c1 = (floor((i + 1 + dimX) / dimX) - 1) * (dimX - 1) +
           (i + 1 + dimX) % dimX;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // symmetry down
  fprintf(file, "(13 (8 %x %x 7 2)(\n",
          numFluentFacesInterior + 2 * numFluentFacesIO + numFluentFacesUp +
              numFluentFacesDiamond + 1,
          numFluentFacesInterior + 2 * numFluentFacesIO + numFluentFacesUp +
              numFluentFacesDiamond + numFluentFacesDown);
  for (int i = 0; i < numFluentNodes; ++i) {
    row = floor(i / dimX);
    col = i % dimX;

    if (row == 0 && col < dimX - 1 &&
        fNodes[i].x > meshConfig->sUp + diamond->l) {
      n1 = i + 1;
      n2 = i + 2;
      c1 = (floor((i + 1 + dimX) / dimX) - 1) * (dimX - 1) +
           (i + 1 + dimX) % dimX;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  fprintf(file, "))\n");
  fprintf(file, "\n");

  // top wall
  fprintf(file, "(13 (9 %x %x 3 2)(\n",
          numFluentFacesInterior + 2 * numFluentFacesIO + numFluentFacesUp +
              numFluentFacesDiamond + numFluentFacesDown + 1,
          numFluentFacesInterior + 2 * numFluentFacesIO + numFluentFacesUp +
              numFluentFacesDiamond + numFluentFacesDown + numFluentFacesTop);
  for (int i = 0; i < numFluentNodes; ++i) {
    row = floor(i / dimX);
    col = i % dimX;

    if (row == dimY - 1 && col < dimX - 1) {
      n1 = i + 2;
      n2 = i + 1;
      c1 = (floor(i / dimX) - 1) * (dimX - 1) + i % dimX;
      c2 = 0;
      fprintf(file, "%x %x %x %x\n", n1, n2, c1, c2);
    }
  }
  fprintf(file, "))\n");
  fprintf(file, "\n");
  fprintf(file, "(10 (1 1 %x 1 2)(\n", numFluentNodes);
  for (int i = 0; i < numFluentNodes; ++i) {
    fprintf(file, "%.6e %.6e\n", fNodes[i].x, fNodes[i].y);
  }
  fprintf(file, "))\n");
  fprintf(file, "\n");

  fclose(file);
  printf("FILE WRITTEN!\n");
  free(fNodes);
  free(meshStructure);
  return 0;
}
