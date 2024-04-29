#ifndef CONVERTER_H
#define CONVERTER_H

typedef struct {
  int tag;
  double x, y, z;
} Point; 

typedef struct {
  int tag;
  double minX, minY, minZ;
  double maxX, maxY, maxZ;
  int tagsBoundingPoints[2]; // lines
} Curve;

typedef struct {
  int tag;
  double minX, minY, minZ;
  double maxX, maxY, maxZ;
  int tagsBoundingCurves[4]; // quads
} Surface;

typedef struct {
  Point* points;
  Curve* curves;
  Surface* surfaces;
} Entities;

void readGmsh(const char* fileName, Entities* entities);



#endif // CONVERTER_H
