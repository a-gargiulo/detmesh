#include "parsing.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "diamond.h"
#include "error.h"
#include "mesh.h"

void trim(char** str, int* err) {

  if (*str == NULL) {
    *err = ERROR_NULL_POINTER;
    return;
  }

  // Leading white space
  while (isspace(**str))
    (*str)++;

  if (**str == '\0') {
    *err = SUCCESS;
    return;
  }

  // Trailing whitespace
  char* end = *str + strlen(*str) - 1;
  while (end > *str && isspace(*end))
    end--;
  
  *(end + 1) = '\0';

  *err = SUCCESS;
  return;

}

void format_var(char* str) {
  int len = strlen(str);
  for (int i = 0; i < len; i++) {
    str[i] = tolower(str[i]);
    if (isspace(str[i])) {
      str[i] = '_';
    }
  }
}


void get_var_val(char** buffer, char** var, char** val, int* err) {
    *var = strtok(*buffer, "(=");
    strtok(NULL, "(=");
    *val = strtok(NULL, "(=");
    trim(var, err);
    trim(val, err);
    format_var(*var);
}


int parse_input(const char* filename, Diamond* diamond, double* x_guess, MeshConfig* config) {

  int err;

  char buffer[BUFFER_SIZE];
  char* pBuffer;

  char* val;
  char* var;

  char* guess_list;
  int list_counter;

  FILE* file;
  file = fopen(filename, "r");

  if (file == NULL) {
    logError("Could not open the input file.", ERROR_COULD_NOT_OPEN_FILE);
    return ERROR_COULD_NOT_OPEN_FILE;
  }


  while (fgets(buffer, BUFFER_SIZE, file) != NULL) {
    pBuffer = buffer;
    trim(&pBuffer, &err);
    if (*pBuffer == '#' || *pBuffer == '\0') {
      continue;
    }

    get_var_val(&pBuffer, &var, &val, &err);

    // Diamond Geometry
    if (strcmp(var, "leading_half_angle") == 0)
      diamond->alpha = atof(val);
    else if (strcmp(var, "trailing_half_angle") == 0)
      diamond->beta = atof(val);
    else if (strcmp(var, "length_of_diamond") == 0)
      diamond->l = strtod(val, NULL);
    else if (strcmp(var, "apex_radius_of_curvature") == 0)
      diamond->r = strtod(val, NULL);
    // Initial Guess
    else if (strcmp(var, "initial_guess") == 0) {
      guess_list = strtok(val, "[,]");
      list_counter = 0;
      while (guess_list != NULL) {
        x_guess[list_counter] = atof(guess_list) / 100.0 * diamond->l;
        guess_list = strtok(NULL, "[,]");
        list_counter++;
      }
    }
    // Mesh configuration
    else if (strcmp(var, "width_upstream_approach_region") == 0)
      config->sUp = strtod(val, NULL);
    else if (strcmp(var, "width_downstream_wake_region") == 0)
      config->sDown = strtod(val, NULL);
    else if (strcmp(var, "width_diamond_leading_edge_cluster") == 0)
      config->sBlkUp = strtod(val, NULL);
    else if (strcmp(var, "width_diamond_trailing_edge_cluster") == 0)
      config->sBlkDown = strtod(val, NULL);
    else if (strcmp(var, "width_upstream_portion_apex_arc_cluster") == 0)
      config->sArcUp = strtod(val, NULL);
    else if (strcmp(var, "width_downstream_portion_apex_arc_cluster") == 0)
      config->sArcDown = strtod(val, NULL);
    else if (strcmp(var, "tunnel_half_height") == 0)
      config->tHeight = strtod(val, NULL);
  }

  fclose(file);
  return 0;
}
