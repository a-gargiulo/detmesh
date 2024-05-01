#include "converter.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_BUFFER_SIZE 1024

void print_title() {
  const char* title = "\n"
                      "READING GMSH\n"
                      "------------\n"
                      "\n";
  printf("%s", title);
}

int readGmsh(const char *fileName, Point **points, int *numPoints,
              Curve **curves, int *numCurves, Surface **surfaces,
              int *numSurfaces, Node **nodes, int *numEntityBlocks,
              Element **elements, int *numEntityBlocksElem, int* numNodes) {

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
  }
  else {
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
        fprintf(stderr, "ERROR: Failed to allocate memory for Points, Curves, or Surfaces!\n");
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

        if ((*nodes)[i].nodeTags == NULL || (*nodes)[i].x == NULL || (*nodes)[i].y == NULL || (*nodes)[i].z == NULL) {
          fprintf(stderr, "ERROR: Failed to allocate memory for node tags or x, y, z coordinates!\n");
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
          fprintf(stderr, "ERROR: Failed to allocate memory for element tags!\n");
          return 1;
        }

        if ((*elements)[i].elementType == 15) {
          (*elements)[i].nodeTags =
              (int *)malloc((*elements)[i].numElementsInBlock * sizeof(int));
          
          if ((*elements)[i].nodeTags == NULL) {
            fprintf(stderr, "ERROR: Failed to allocate memory for node tags!\n");
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
            fprintf(stderr, "ERROR: Failed to allocate memory for node tags!\n");
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
            fprintf(stderr, "ERROR: Failed to allocate memory for node tags!\n");
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
  else
  {
    fprintf(stderr, "ERROR: Could not close file!\n");
    return 1; 
  }
  return 0;
}


int writeFluent(const char* outputFile, const Node* nodes, const int numEntityBlocks, const int numNodes) {
  FILE* file;


  double* sNodeX = (double*)malloc(numNodes * sizeof(double));  
  double* sNodeY = (double*)malloc(numNodes * sizeof(double));  
  int* sNodeTags = (int*)malloc(numNodes * sizeof(int));  
  if (sNodeX == NULL || sNodeY == NULL) {
    fprintf(stderr, "Error: Memory allocation failed\n");
    // Handle memory allocation failure
  }

  // flatten nodes
  int k = 0;
  for (int i = 0; i < numEntityBlocks; ++i) {
    for (int j = 0; j < nodes[i].numNodesInBlock; ++j)
    {
      sNodeX[k] = nodes[i].x[j];
      sNodeY[k] = nodes[i].y[j];
      sNodeTags[k] = nodes[i].nodeTags[j];
      k++;
    }
  }

  // sort nodes



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
  fprintf(file, "(10 (0 1 %x 0 2))\n", numNodes);
  fprintf(file, "\n");
  fprintf(file, "(10 (1 1 %x 1 2)(\n", numNodes);

  fprintf(file, ")\n");

  fclose(file);
  free(sNodeX);
  free(sNodeY);
}


