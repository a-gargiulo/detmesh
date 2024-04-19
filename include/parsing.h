#ifndef PARSING_H
#define PARSING_H

#define BUFFER_SIZE 100

#include "diamond.h"

int trim(char** str);

int firstNonSpaceIndex(const char* str);

int prepare_var(char* str);

int read_input(const char* filename, Diamond* diamond, double* x_guess);

#endif  // PARSING_H
