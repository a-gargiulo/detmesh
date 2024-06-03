#ifndef PARSING_H
#define PARSING_H

#include <stdbool.h>
#include <stddef.h>

#include "diamond.h"
#include "gmesh.h"


typedef struct {
  const char* name;
  double* ptr_to_value;
} VariableMapping;

void trim(char** str);

void format_variable_name(char* str);

bool is_variable(const char* var, const char* name);

int update_variable(const char* name, const double* value, VariableMapping* mapping, size_t n_mapping);

int parse_user_input(const char* file_name, Diamond* diamond, double* x_guess,
                     GMeshConfig* config);

#endif  // PARSING_H
