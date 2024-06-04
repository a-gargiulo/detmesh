#include "error.h"

#include <stdio.h>

void log_error(const char* message, int error_code)
{
    fprintf(stderr, "Error %d: %s\n", error_code, message);
}
