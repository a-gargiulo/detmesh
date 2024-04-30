#include "converter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_BUFFER_SIZE 1024

void readGmsh(const char *fileName, Point **points, size_t *numPoints,
              Curve **curves, size_t *numCurves, Surface **surfaces,
              size_t *numSurfaces, Node **nodes, size_t *numEntityBlocks,
              Element **elements, size_t *numEntityBlocksElem) {
  int err;
  int isFirstLine;
  int readPts, readCrvs, readSurfs;
  int iiPts, iiCrvs, iiSurfs;
  size_t numNodes, minNodeTag, maxNodeTag;
  size_t numElements, minElementTag, maxElementTag;
  FILE *file;
  char line[LINE_BUFFER_SIZE];
  char *pLine;
  char startCategory[50];
  char endCategory[50];

  printf("\n");

  printf("Converter\n");

  printf("---------\n");

  printf("\n");

  file = fopen(fileName, "r");

  if (file == NULL) {
    fprintf(stderr, "ERROR: Gmsh file could not open!\n");

    return;
  }

  while (fgets(line, LINE_BUFFER_SIZE, file) != NULL || !feof(file)) {
    while (sscanf(line, " $%s ", startCategory) != 1) {
      fgets(line, LINE_BUFFER_SIZE, file);
    }
    printf("%s\n", startCategory);

    if (strcmp(startCategory, "Entities") == 0) {
      printf("Reading ENTITIES... ");

      // first line
      fgets(line, LINE_BUFFER_SIZE, file);
      sscanf(line, " %zu %zu %zu 0", numPoints, numCurves, numSurfaces);
      *points = (Point *)malloc(*numPoints * sizeof(Point));
      *curves = (Curve *)malloc(*numCurves * sizeof(Curve));
      *surfaces = (Surface *)malloc(*numSurfaces * sizeof(Curve));

      // points
      for (size_t i = 0; i < *numPoints; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %lf %lf %lf ", &(*points)[i].tag, &(*points)[i].x,
               &(*points)[i].y, &(*points)[i].z);
      }

      // curves
      for (size_t i = 0; i < *numCurves; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %lf %lf %lf %lf %lf %lf 0 2 %d %d", &(*curves)[i].tag,
               &(*curves)[i].minX, &(*curves)[i].minY, &(*curves)[i].minZ,
               &(*curves)[i].maxX, &(*curves)[i].maxY, &(*curves)[i].maxZ,
               &(*curves)[i].tagsBoundingPoints[0],
               &(*curves)[i].tagsBoundingPoints[1]);
      }

      // surfaces
      for (size_t i = 0; i < *numSurfaces; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %lf %lf %lf %lf %lf %lf 0 4 %d %d %d %d",
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
        printf("DONE!\n");
      }
    } else if (strcmp(startCategory, "Nodes") == 0) {
      printf("Reading NODES... ");

      // first line
      fgets(line, LINE_BUFFER_SIZE, file);
      sscanf(line, "%zu %zu %zu %zu", numEntityBlocks, &numNodes, &minNodeTag,
             &maxNodeTag);
      *nodes = (Node *)malloc(*numEntityBlocks * sizeof(Node));

      // entity blocks
      for (size_t i = 0; i < *numEntityBlocks; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %d 0 %zu ", &(*nodes)[i].entityDim,
               &(*nodes)[i].entityTag, &(*nodes)[i].numNodesInBlock);
        (*nodes)[i].nodeTags =
            (size_t *)malloc((*nodes)[i].numNodesInBlock * sizeof(size_t));
        (*nodes)[i].x =
            (double *)malloc((*nodes)[i].numNodesInBlock * sizeof(double));
        (*nodes)[i].y =
            (double *)malloc((*nodes)[i].numNodesInBlock * sizeof(double));
        (*nodes)[i].z =
            (double *)malloc((*nodes)[i].numNodesInBlock * sizeof(double));

        for (size_t j = 0; j < (*nodes)[i].numNodesInBlock; ++j) {
          fgets(line, LINE_BUFFER_SIZE, file);
          sscanf(line, " %zu ", &(*nodes)[i].nodeTags[j]);
        }

        for (size_t k = 0; k < (*nodes)[i].numNodesInBlock; ++k) {
          fgets(line, LINE_BUFFER_SIZE, file);
          sscanf(line, " %lf %lf %lf ", &(*nodes)[i].x[k], &(*nodes)[i].y[k],
                 &(*nodes)[i].z[k]);
        }
      }

      while (sscanf(line, " $%s ", endCategory) != 1) {
        fgets(line, LINE_BUFFER_SIZE, file);
      }

      if (strcmp(endCategory, "EndNodes") == 0) {
        printf("DONE!\n");
      }
    } else if (strcmp(startCategory, "Elements") == 0) {
      printf("Reading ELEMENTS... ");

      // first line
      fgets(line, LINE_BUFFER_SIZE, file);
      sscanf(line, "%zu %zu %zu %zu", numEntityBlocksElem, &numElements,
             &minElementTag, &maxElementTag);
      *elements = (Element *)malloc(*numEntityBlocksElem * sizeof(Element));

      // Element blocks
      for (size_t i = 0; i < *numEntityBlocksElem; ++i) {
        fgets(line, LINE_BUFFER_SIZE, file);
        sscanf(line, " %d %d %d %zu ", &(*elements)[i].entityDim,
               &(*elements)[i].entityTag, &(*elements)[i].elementType,
               &(*elements)[i].numElementsInBlock);
        (*elements)[i].elementTags = (size_t *)malloc(
            (*elements)[i].numElementsInBlock * sizeof(size_t));
        if ((*elements)[i].elementType == 15) {
          (*elements)[i].nodeTags = (size_t *)malloc(
              (*elements)[i].numElementsInBlock * sizeof(size_t));
          for (size_t j = 0; j < (*elements)[i].numElementsInBlock; ++j) {
            fgets(line, LINE_BUFFER_SIZE, file);
            sscanf(line, " %zu %zu ", &(*elements)[i].elementTags[j],
                   &(*elements)[i].nodeTags[j]);
          }
        } else if ((*elements)[i].elementType == 1) {
          (*elements)[i].nodeTags = (size_t *)malloc(
              2 * (*elements)[i].numElementsInBlock * sizeof(size_t));
          for (size_t j = 0; j < (*elements)[i].numElementsInBlock; ++j) {
            fgets(line, LINE_BUFFER_SIZE, file);
            sscanf(line, " %zu %zu %zu ", &(*elements)[i].elementTags[j],
                   &(*elements)[i].nodeTags[j * 2],
                   &(*elements)[i].nodeTags[j * 2 + 1]);
          }
        } else if ((*elements)[i].elementType == 3) {
          (*elements)[i].nodeTags = (size_t *)malloc(
              4 * (*elements)[i].numElementsInBlock * sizeof(size_t));
          for (size_t j = 0; j < (*elements)[i].numElementsInBlock; ++j) {
            fgets(line, LINE_BUFFER_SIZE, file);
            sscanf(line, "%zu %zu %zu %zu %zu ", &(*elements)[i].elementTags[j],
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
        printf("DONE!\n");
      }
    }
  }
  if (feof(file)) {
    printf("END OF FILE!\n");
  }
  fclose(file);
}
