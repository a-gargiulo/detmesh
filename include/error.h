#ifndef ERROR_H
#define ERROR_H

#define SUCCESS 0
#define ERROR_NULL_POINTER -1
#define ERROR_INVALID_INPUT -2
#define ERROR_OUT_OF_MEMORY -3
#define ERROR_COULD_NOT_OPEN_FILE -4
#define ERROR_VARIABLE_NOT_FOUND -5
#define ERROR_DIVISION_BY_ZERO -6

void log_error(const char *message, int error_code);

#endif // ERROR_H
