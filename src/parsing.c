#include "parsing.h"

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diamond.h"
#include "error.h"
#include "gmesh.h"

#define PARSER_LINE_BUFFER_SIZE 1024
#define PARSER_VARIABLE_NAME_BUFFER_SIZE 100

void trim(char **str)
{
    while (isspace(**str))
        (*str)++;

    if (**str == '\0')
    {
        printf("WARNING: The input string is empty!\n");
    }

    char *end = *str + strlen(*str) - 1;
    while (end > *str && isspace(*end))
        end--;

    *(end + 1) = '\0';
}

void format_variable_name(char *str)
{
    int len = strlen(str);
    for (int i = 0; i < len; i++)
    {
        str[i] = tolower(str[i]);
        if (isspace(str[i]))
        {
            str[i] = '_';
        }
    }
}

bool is_variable(const char *var, const char *name)
{
    return strcmp(var, name) == 0;
}

int update_variable(const char *name, const double *value, VariableMapping *mapping, size_t n_mapping)
{
    for (size_t i = 0; i < n_mapping; ++i)
    {
        if (is_variable(name, mapping[i].name))
        {
            *(mapping[i].ptr_to_value) = *value;
            return SUCCESS;
        }
    }
    log_error("Variable not found!", ERROR_VARIABLE_NOT_FOUND);
    return ERROR_VARIABLE_NOT_FOUND;
}

int parse_user_input(const char *file_name, Diamond *diamond, double *x_guess, GMeshConfig *mesh_config)
{
    int status;

    VariableMapping mapping[] = {{"leading_half_angle", &diamond->alpha},
                                 {"trailing_half_angle", &diamond->beta},
                                 {"length_of_diamond", &diamond->l},
                                 {"apex_radius_of_curvature", &diamond->r},
                                 {"width_upstream_approach_region", &mesh_config->sUp},
                                 {"width_downstream_wake_region", &mesh_config->sDown},
                                 {"width_diamond_leading_edge_cluster", &mesh_config->sBlkUp},
                                 {"width_diamond_trailing_edge_cluster", &mesh_config->sBlkDown},
                                 {"width_upstream_portion_apex_arc_cluster", &mesh_config->sArcUp},
                                 {"width_downstream_portion_apex_arc_cluster", &mesh_config->sArcDown},
                                 {"tunnel_half_height", &mesh_config->tHeight}};

    char line_buffer[PARSER_LINE_BUFFER_SIZE];
    char name_buffer[PARSER_VARIABLE_NAME_BUFFER_SIZE];
    char *name;
    double value;

    FILE *file = fopen(file_name, "r");
    if (file == NULL)
    {
        log_error("Could not open the input file.", ERROR_COULD_NOT_OPEN_FILE);
        return ERROR_COULD_NOT_OPEN_FILE;
    }

    while (fgets(line_buffer, PARSER_LINE_BUFFER_SIZE, file) != NULL && !feof(file))
    {
        if (sscanf(line_buffer, " %[^(] (%*[^)]) = %lf ", name_buffer, &value) == 2)
        {
            name = name_buffer;
            trim(&name);
            format_variable_name(name);
            status = update_variable(name, &value, mapping, sizeof(mapping) / sizeof(mapping[0]));
            if (status != 0)
            {
                return status;
            }
        }
        else if (sscanf(line_buffer, " %[^(] (%*[^)]) = [ %lf , %lf , %lf , %lf ]", name_buffer, &x_guess[0], &x_guess[1],
                        &x_guess[2], &x_guess[3]) == 5)
        {
            for (size_t i = 0; i < 4; ++i)
            {
                x_guess[i] = x_guess[i] / 100.0 * diamond->l;
            }
        }
    }

    fclose(file);
    return 0;
}
