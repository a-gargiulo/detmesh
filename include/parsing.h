#ifndef PARSING_H
#define PARSING_H

#define BUFFER_SIZE 100

#include "diamond.h"
#include "mesh.h"

void trim(char** str, int* err);

void get_var_val(char** buffer, char** var, char** val, int* err);

void format_var(char* str);

int parse_input(const char* filename, Diamond* diamond, double* x_guess, MeshConfig* config);

#endif  // PARSING_H
