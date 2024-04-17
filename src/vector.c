#include "utility.h"

#include <math.h>

Vector* create_vector(int size) {
  Vector* vec = (Vector*)malloc(sizeof(Vector));
  if (vec == NULL) {
    fprintf()
  }
  vec->size = size;
  vec->data = (double*)malloc(size * sizeof(double));

  return vec;
}

double vec_norm(const double* vec, int size) {
  double sum = 0.0;
  for (int i = 0; i < size; ++i) {
    sum += vec[i] * vec[i];
  }
  return sqrt(sum);
}

void vec_subtract(const double* vec1, const double* vec2, double* result,
                  int size) {
  for (int i = 0; i < size; ++i) {
    result[i] = vec1[i] - vec2[i];
  }
}

void vec_division(const double* vec1, const double* vec2, double* result,
                  int size) {
  for (int i = 0; i < size; ++i) {
    result[i] = vec1[i] / vec2[i];
  }
}
