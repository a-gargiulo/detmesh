#include "converter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include "parsing.h"

#define LINE_BUFFER_SIZE 1024

void readGmsh(const char* fileName, Point** points, size_t* numPoints, Curve** curves, size_t* numCurves, Surface** surfaces, size_t* numSurfaces, Node** nodes) {
  int err;
  int isFirstLine;
  int readPts, readCrvs, readSurfs, readVols;
  int iiPts, iiCrvs, iiSurfs;
  size_t numEntityBlocks, numNodes, minNodeTag, maxNodeTag; 
  FILE* file;
  char line[LINE_BUFFER_SIZE];
  char* pLine;
  char startCategory[50];
  char endCategory[50];

  file = fopen(fileName, "r");

  if (file == NULL) {
    fprintf(stderr, "MESH FILE COULD NOT OPEN\n");
    return;
  }

  while (fgets(line, LINE_BUFFER_SIZE, file) != NULL || !feof(file)) {
    pLine = line;
    trim(&pLine, &err);

    if (sscanf(line, " $%s", startCategory) == 1) {
      sprintf(endCategory, "$End%s", startCategory);
      isFirstLine = 1;
      continue;
    }

    if (strcmp(pLine, endCategory) != 0) {

      if (strcmp(startCategory, "Entities") == 0) {
        if (isFirstLine) {
          sscanf(line, " %zu %zu %zu 0", numPoints, numCurves, numSurfaces); 
          *points = (Point*)malloc(*numPoints * sizeof(Point));
          *curves = (Curve*)malloc(*numCurves * sizeof(Curve));
          *surfaces = (Surface*)malloc(*numSurfaces * sizeof(Curve));
          iiPts = 0;
          iiCrvs = 0;
          iiSurfs = 0;
          readPts = 1;
          isFirstLine = 0;
        } 
        else {
          if (readPts!=0) {
            sscanf(line, " %d %lf %lf %lf ", &(*points)[iiPts].tag, &(*points)[iiPts].x, &(*points)[iiPts].y, &(*points)[iiPts].z);
            iiPts++;
            if (iiPts == *numPoints) {
              readPts = 0;
              readCrvs = 1;
            }
          }
          else if (readCrvs!=0){
            sscanf(line, " %d %lf %lf %lf %lf %lf %lf 0 2 %d %d", &(*curves)[iiCrvs].tag, &(*curves)[iiCrvs].minX, &(*curves)[iiCrvs].minY, &(*curves)[iiCrvs].minZ, &(*curves)[iiCrvs].maxX, &(*curves)[iiCrvs].maxY, &(*curves)[iiCrvs].maxZ, &(*curves)[iiCrvs].tagsBoundingPoints[0], &(*curves)[iiCrvs].tagsBoundingPoints[1]);
            iiCrvs++;
            if (iiCrvs == *numCurves) {
              readCrvs = 0;
              readSurfs = 1;
            }
          }
          else if ( readSurfs!=0 ) {
            sscanf(line, " %d %lf %lf %lf %lf %lf %lf 0 4 %d %d %d %d", &(*surfaces)[iiSurfs].tag, &(*surfaces)[iiSurfs].minX, &(*surfaces)[iiSurfs].minY, &(*surfaces)[iiSurfs].minZ, &(*surfaces)[iiSurfs].maxX, &(*surfaces)[iiSurfs].maxY, &(*surfaces)[iiSurfs].maxZ, &(*surfaces)[iiSurfs].tagsBoundingCurves[0], &(*surfaces)[iiSurfs].tagsBoundingCurves[1], &(*surfaces)[iiSurfs].tagsBoundingCurves[2], &(*surfaces)[iiSurfs].tagsBoundingCurves[3]);
            iiSurfs++;
            if (iiSurfs == *numSurfaces) {
              readSurfs = 0;
            }
          }
        }
      }

      // else if (strcmp(startCategory, "Nodes") == 0) {

      // }




    }
  }

  // printf("STATUS %d\n", stat);
  // getchar();
  fclose(file);
  // fclose(file);
}
