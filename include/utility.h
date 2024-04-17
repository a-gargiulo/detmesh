#ifndef UTILITY_H
#define UTILITY_H

void vec_norm(const double* vec, double* result, int size);

// Vector subtraction: result = vec1 - vec2
void vec_subtract(const double* vec1, const double* vec2, double* result,
                  int size);

// Element-wise division: result = vec1 / vec2
void vec_division(const double* vec1, const double* vec2, double* result,
                  int size);

#endif  // UTILITY_H
