#include "parsing.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diamond.h"

char* trim(const char* str) {
  int original_length = strlen(str);

  int start = 0;
  while (str[start] == ' ') {
    start++;
  }

  int end = original_length - 1;
  while (str[end] == ' ') {
    end--;
  }

  int new_length = end - start + 1;
  char* trimmed_str = (char*)malloc((new_length + 1) * sizeof(char));

  int j = 0;
  for (int i = start; i <= end; ++i) {
    trimmed_str[j] = str[i];
    j++;
  }
  trimmed_str[j] = '\0';
  return trimmed_str;
}

char* prepare_var(const char* str) {
  int len = strlen(str);
  char* lower_str = (char*)malloc((len + 1) * sizeof(char));
  strcpy(lower_str, str);
  for (int i = 0; i < len; i++) {
    lower_str[i] = tolower(lower_str[i]);
    if (lower_str[i] == ' ') {
      lower_str[i] = '_';
    }
  }

  return lower_str;
}

int read_input(const char* filename, Geometry* geo, double* x) {
  FILE* pFile;
  char line[100];
  pFile = fopen(filename, "r");
  if (pFile == NULL) {
    fprintf(stderr, "ERROR: Could not open file!\n");
    return 1;
  }

  char* var;
  char* num;
  char* token;
  char* delim = "=\n";
  while (fgets(line, sizeof(line), pFile) != NULL) {
    char* trim_line = trim(line);
    if (trim_line[0] == '#' || trim_line[0] == '\n') {
      continue;
    }
    token = strtok(line, delim);
    num = strtok(NULL, delim);
    var = strtok(token, "(");

    char* trim_var = trim(var);
    char* trim_num = trim(num);
    char* prep_var = prepare_var(trim_var);

    // printf("%s\n", prep_var);
    // printf("%s\n", trim_num);
    // printf("%d\n", strcmp(prep_var, "trailing_half_angle"));
    if (strcmp(prep_var, "leading_half_angle") == 0) {
      geo->alpha = atof(trim_num);
    } else if (strcmp(prep_var, "trailing_half_angle") == 0) {
      geo->beta = atof(trim_num);
    } else if (strcmp(prep_var, "apex_radius_of_curvature") == 0) {
      char* endptr1;
      geo->r = strtod(trim_num, &endptr1);
    } else if (strcmp(prep_var, "length_of_diamond") == 0) {
      char* endptr2;
      geo->l_d = strtod(trim_num, &endptr2);
    } else if (strcmp(prep_var, "initial_guess") == 0) {
      char* list;
      char* delims = "[,]";
      list = strtok(trim_num, delims);
      int counter = 0;
      while (list != NULL) {
        x[counter] = atof(list) / 100.0 * geo->l_d;
        list = strtok(NULL, "[,]");
        counter++;
      }
    }
    free(trim_var);
    free(trim_num);
    free(trim_line);
    free(prep_var);  // }
  }
  fclose(pFile);
  return 0;
}
