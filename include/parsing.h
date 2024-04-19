#ifndef PARSING_H
#define PARSING_H

#define BUFFER_SIZE 100

#include "diamond.h"

int trim(char** str);

int get_first_non_space_index(const char* str);

int format_variable_name(char* str);

int read_input(const char* filename, Diamond* diamond, double* x_guess);

#endif  // PARSING_H
