#include "converter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parsing.h"

#define LINE_BUFFER_SIZE 1024

void readGmsh(const char* fileName, Entities* entities) {
  int err;
  int isFirstLine;
  int readPts, readCrvs, readSurfs, readVols;
  int iiPts, iiCrvs, iiSurfs;
  size_t numPoints, numCurves, numSurfaces, numVolumes; 
  FILE* file;
  char line[LINE_BUFFER_SIZE];
  char* pLine;
  char startCategory[50];
  char endCategory[50];

  file = fopen(fileName, "r");

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
          sscanf(line, " %zu %zu %zu %zu ", &numPoints, &numCurves, &numSurfaces, &numVolumes); 
          Point* points = (Point*)malloc(numPoints * sizeof(Point));
          Curve* curves = (Curve*)malloc(numCurves * sizeof(Curve));
          Surface* surfaces = (Surface*)malloc(numSurfaces * sizeof(Curve));
          entities->points = points;
          entities->curves = curves;
          entities->surfaces = surfaces;  
          iiPts = 0;
          iiCrvs = 0;
          iiSurfs = 0;
          readPts = 1;
          isFirstLine = 0;
        } 
        else {
          if (readPts!=0) {
            sscanf(line, " %d %lf %lf %lf ", &(entities->points[iiPts].tag), &(entities->points[iiPts].x), &(entities->points[iiPts].y), &(entities->points[iiPts].z));
            iiPts++;
            if (iiPts == numPoints) {
              readPts = 0;
              readCrvs = 1;
            }
          }
          else if (readCrvs!=0){
            sscanf(line, " %d %lf %lf %lf %lf %lf %lf 0 2 %d %d", &(entities->curves[iiCrvs].tag), &(entities->curves[iiCrvs].minX), &(entities->curves[iiCrvs].minY), &(entities->curves[iiCrvs].minZ), &(entities->curves[iiCrvs].maxX), &(entities->curves[iiCrvs].maxY), &(entities->curves[iiCrvs].maxZ), &(entities->curves[iiCrvs].tagsBoundingPoints[0]), &(entities->curves[iiCrvs].tagsBoundingPoints[1]));
            iiCrvs++;
            if (iiCrvs == numCurves) {
              readCrvs = 0;
              readSurfs = 1;
            }
          }
          else if ( readSurfs!=0 ) {
            sscanf(line, " %d %lf %lf %lf %lf %lf %lf 0 4 %d %d %d %d", &(entities->surfaces[iiSurfs].tag), &(entities->surfaces[iiSurfs].minX), &(entities->surfaces[iiSurfs].minY), &(entities->surfaces[iiSurfs].minZ), &(entities->surfaces[iiSurfs].maxX), &(entities->surfaces[iiSurfs].maxY), &(entities->surfaces[iiSurfs].maxZ), &(entities->surfaces[iiSurfs].tagsBoundingCurves[0]), &(entities->surfaces[iiSurfs].tagsBoundingCurves[1]), &(entities->surfaces[iiSurfs].tagsBoundingCurves[2]), &(entities->surfaces[iiSurfs].tagsBoundingCurves[3]));
            iiSurfs++;
            if (iiSurfs == numSurfaces) {
              readSurfs = 0;
            }
          }
        }
      }
    }
  }

  fclose(file);

  for (size_t i = 0; i < numPoints; ++i) {
    printf("%d %.12lf %.12lf %.12lf\n", (entities->points[i].tag), (entities->points[i].x), (entities->points[i].y), (entities->points[i].z));
  }
}
