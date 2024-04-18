#ifndef PARSING_H
#define PARSING_H

#include "diamond.h"

char* trim(const char* str);
char* prepare_var(const char* str);
int read_input(const char* filename, Geometry* geo, double* x);

#endif  // PARSING_H
