#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vector.h"

Vector* construct_vector(const int size) {
  Vector* vec = (Vector*)malloc(sizeof(Vector));
  if (vec == NULL) {
    fprintf(stderr, "ERROR! Could not allocate memory for vector.\n");
    return NULL;
  }
  vec->size = size;
  vec->data = (double*)calloc(size, sizeof(double));
  if (vec->data == NULL) {
    free(vec);
    fprintf(stderr, "ERROR! Could not allocate memory for vector data.\n");
    return NULL;
  }
  return vec;
}

int delete_vector(Vector* vec) {
  if (vec == NULL) {
    fprintf(stderr, "ERROR! Attempting to delete NULL vector.\n");
    return -1;
  }

  if (vec->data != NULL) {
    memset(vec->data, 0, vec->size * sizeof(double));
  }
  free(vec->data);
  free(vec);
  return 0;
}

double vector_norm(const Vector* vec) {
  double sum = 0.0;
  for (int i = 0; i < vec->size; ++i) {
    sum += vec->data[i] * vec->data[i];
  }
  return sqrt(sum);
}

void vector_subtract(const Vector* vec1, const Vector* vec2, Vector* result) {
  if (vec1->size != vec2->size || vec1->size != result->size ||
      vec2->size != result->size) {
    fprintf(stderr, "Error! Vectors must be of the same size.\n");
    // return NULL;
  }
  for (int i = 0; i < result->size; ++i) {
    result->data[i] = vec1->data[i] - vec2->data[i];
  }
}

void vector_divide(const Vector* vec1, const Vector* vec2, Vector* result) {
  if (vec1->size != vec2->size || vec1->size != result->size ||
      vec2->size != result->size) {
    fprintf(stderr, "Error: Vectors must be of the same size.\n");
    // return NULL;
  }
  for (int i = 0; i < result->size; ++i) {
    result->data[i] = vec1->data[i] / vec2->data[i];
  }
}
