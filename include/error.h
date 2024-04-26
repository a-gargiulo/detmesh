#ifndef ERROR_H
#define ERROR_H

#define SUCCESS 0
#define ERROR_NULL_POINTER -1
#define ERROR_INVALID_INPUT -2
#define ERROR_OUT_OF_MEMORY -3
#define ERROR_COULD_NOT_OPEN_FILE -4

void logError(const char *message, int errorCode);

#endif // ERROR_H