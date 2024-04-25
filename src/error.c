#include "error.h"

#include <stdio.h>

void logError(const char *message, int errorCode) {
    fprintf(stderr, "Error %d: %s\n", errorCode, message);
}
