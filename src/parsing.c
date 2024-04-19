#include "parsing.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diamond.h"

int trim(char** str) {
  // Works with string literals
  if (*str == NULL) {
    printf("Error: String pointer is NULL.\n");
    return 1;
  }

  // Trim leading whitespace
  while (isspace(**str)) {
    (*str)++;
  }

  if (**str == '\0') {
    **str = '\0';
    return 0;
  }

  // Trim trailing whitespace
  char *end = *str + strlen(*str) - 1;
  while (end > *str && isspace(*end)) {
  // while (end > *str && *end == '-') {
    end--;
  }
  *(end+1) = '\0';

  return 0;
}


int prepare_var(char* str) {
  int len = strlen(str);
  for (int i = 0; i < len; i++) {
    str[i] = tolower(str[i]);
    if (isspace(str[i])) {
      str[i] = '_';
    }
  }

  return 0;
}


int firstNonSpaceIndex(const char *str) {
    int index = 0;
    if (str[index] == '\n')
    {
      return index;
    }
    while (str[index] != '\0' && isspace(str[index])) {
        index++;
    }
    return index;
}


int read_input(const char* filename, Geometry* geo, double* x_guess) {
  FILE* pFile;
  pFile = fopen(filename, "r");

  if (pFile == NULL) {
    fprintf(stderr, "ERROR: Could not open file!\n");
    return 1;
  }

  char buffer[BUFFER_SIZE];

  char* tmp;
  char* num;
  char* var;
  while (fgets(buffer, BUFFER_SIZE, pFile) != NULL) {
    // trim(buffer);

    int idx = firstNonSpaceIndex(buffer);
    // if (buffer[idx] == '#' || buffer[0]=='\n') {
    if (buffer[idx] == '#'|| buffer[idx]=='\n') {
      continue;
    }

    var = strtok(buffer, "(=");
    strtok(NULL, "(=");
    num = strtok(NULL, "(=");
    trim(&var);
    trim(&num);
    prepare_var(var);

    printf("%s\n", var);
    printf("%s\n", num);
    if (strcmp(var, "leading_half_angle") == 0) {
      geo->alpha = atof(num);
    } else if (strcmp(var, "trailing_half_angle") == 0) {
      geo->beta = atof(num);
    } else if (strcmp(var, "apex_radius_of_curvature") == 0) {
      char* endptr1;
      geo->r = strtod(num, &endptr1);
    } else if (strcmp(var, "length_of_diamond") == 0) {
      char* endptr2;
      geo->l_d = strtod(num, &endptr2);
    } else if (strcmp(var, "initial_guess") == 0) {
      char* list;
      char* delims = "[,]";
      list = strtok(num, delims);
      int counter = 0;
      while (list != NULL) {
        x_guess[counter] = atof(list) / 100.0 * geo->l_d;
        list = strtok(NULL, "[,]");
        counter++;
      }
    }
  }
  fclose(pFile);
  return 0;
}
