#ifndef VECTOR_H
#define VECTOR_H

typedef struct {
  int size;
  double* data;
} Vector;

Vector* construct_vector(const int size);

int delete_vector(Vector* vec);

double vector_norm(const Vector* vec);

void vector_subtract(const Vector* vec1, const Vector* vec2, Vector* result);

void vector_divide(const Vector* vec1, const Vector* vec2, Vector* result);

#endif  // VECTOR_H
