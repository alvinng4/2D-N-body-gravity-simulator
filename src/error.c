#include <stdio.h>

#include "error.h"
#include "gravity_sim.h"

IN_FILE int get_error_msg(const int error_code, char **error_msg);

WIN32DLL_API void print_error_msg(const int error_code)
{
    char *error_msg = NULL;
    int return_code = get_error_msg(error_code, &error_msg);
    if (return_code == SUCCESS)
    {
        fprintf(stderr, "%s", error_msg);
    }
    else
    {
        fprintf(stderr, "Error: error detected but error code is not recognized. (Error Code: %d)\n", error_code);
    }
}

IN_FILE int get_error_msg(const int error_code, char **error_msg)
{
    switch (error_code)
    {
        /* Acceleration */
        // get_acceleration_method_flag()
        case ERROR_UNKNOWN_ACCELERATION_METHOD:
            *error_msg = "Error: Acceleration method not recognized in get_acceleration_method_flag().\n";
            return SUCCESS;

        /* Solution storing (general) */
        // case ERROR_UNKNOWN_STORING_METHOD:
        
        //     break
        default:
            return ERROR_UNKNOWN_ERROR_CODE;
    }
}
