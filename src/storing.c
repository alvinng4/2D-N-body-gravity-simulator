#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "gravity_sim.h"
#include "storing.h"
#include "utils.h"

WIN32DLL_API int get_storing_method_flag(
    const char *restrict storing_method,
    uint *restrict storing_method_flag
)
{
    if (strcmp(storing_method, "default") == 0)
    {
        *storing_method_flag = STORING_METHOD_DEFAULT;
        return SUCCESS;
    }
    else if (strcmp(storing_method, "flush") == 0)
    {
        *storing_method_flag = STORING_METHOD_FLUSH;
        return SUCCESS;
    }
    else if (strcmp(storing_method, "disabled") == 0)
    {
        *storing_method_flag = STORING_METHOD_DISABLED;
        return SUCCESS;
    }
    else
    {
        return ERROR_UNKNOWN_STORING_METHOD;
    }
}

WIN32DLL_API int allocate_solutions_memory(
    Solutions *restrict solutions,
    StoringParam *restrict storing_param,
    int objects_count
)
{
    solutions->sol_state = NULL;
    solutions->sol_time = NULL;
    solutions->sol_dt = NULL;
    
    solutions->sol_state = malloc(storing_param->max_sol_size_ * objects_count * 6 * sizeof(double));
    solutions->sol_time = malloc(storing_param->max_sol_size_ * sizeof(double));
    solutions->sol_dt = malloc(storing_param->max_sol_size_ * sizeof(double));
    if (!solutions->sol_state || !solutions->sol_time || !solutions->sol_dt)
    {
        return ERROR_SOL_OUTPUT_MEMORY_ALLOC;
    }

    return SUCCESS;
}

WIN32DLL_API int open_flush_file(
    StoringParam *restrict storing_param
)
{
    storing_param->flush_file_ = fopen(storing_param->flush_path, "a");
    if (!(storing_param->flush_file_))
    {
        return ERROR_FLUSH_FILE_OPEN;
    }
    return SUCCESS;
}

WIN32DLL_API int close_flush_file(
    StoringParam *restrict storing_param
)
{
    if (storing_param->storing_method_flag_ != STORING_METHOD_FLUSH)
    {
        return ERROR_FLUSH_FILE_CLOSE_NOT_FLUSH_METHOD;
    }
    
    if (!storing_param->flush_file_)
    {
        return ERROR_FLUSH_FILE_CLOSE_IS_NULL;
    }

    if (fclose(storing_param->flush_file_) != 0)
    {
        return ERROR_FLUSH_FILE_CLOSE;
    }
    return SUCCESS;
}

WIN32DLL_API int flush_solution_step_to_csv_file(
    FILE *restrict file,
    const System *restrict system,
    const SimulationStatus *restrict simulation_status,
    Solutions *restrict solutions
)
{
    double energy;
    int return_code;
    int objects_count = system->objects_count;
    real *restrict x = system->x;
    real *restrict v = system->v;

    return_code = compute_energy_step(system, &energy);
    if (return_code != SUCCESS)
    {
        goto error;
    }

    fprintf(file, "%.17g", *(simulation_status->t));
    fprintf(file, ",%.17g", simulation_status->dt);
    fprintf(file, ",%.17g", energy);
    for (int i = 0; i < objects_count; i++)
    {
        fprintf(file, ",%.17g", x[i * 3 + 0]);
        fprintf(file, ",%.17g", x[i * 3 + 1]);
        fprintf(file, ",%.17g", x[i * 3 + 2]);
    }
    for (int i = 0; i < objects_count; i++)
    {
        fprintf(file, ",%.17g", v[i * 3 + 0]);
        fprintf(file, ",%.17g", v[i * 3 + 1]);
        fprintf(file, ",%.17g", v[i * 3 + 2]);
    }
    fprintf(file, "\n");
    fflush(file);

    /* Update solution size */
    *(solutions->sol_size_) += 1;

    return SUCCESS;

error:
    return return_code;
}

WIN32DLL_API int store_solution_step_to_memory(
    const System *restrict system,
    const SimulationStatus *restrict simulation_status,
    Solutions *restrict solutions
)
{
    const int objects_count = system->objects_count;
    real *restrict x = system->x;
    real *restrict v = system->v;
    int64 sol_size = *(solutions->sol_size_);

    /* Store solution */
    memcpy(
        &solutions->sol_state[sol_size * objects_count * 6],
        x,
        objects_count * 3 * sizeof(double)
    );
    memcpy(
        &solutions->sol_state[sol_size * objects_count * 6 + objects_count * 3],
        v,
        objects_count * 3 * sizeof(double)
    );
    solutions->sol_time[sol_size] = *(simulation_status->t);
    solutions->sol_dt[sol_size] = simulation_status->dt;

    /* Update solution size */
    *(solutions->sol_size_) += 1;

    return SUCCESS;
}

WIN32DLL_API int store_solution_step(
    StoringParam *restrict storing_param,
    const System *restrict system,
    const SimulationStatus *restrict simulation_status,
    Solutions *restrict solutions
)
{
    int return_code;
    if (storing_param->storing_method_flag_ == STORING_METHOD_DEFAULT)
    {
        if ((*(solutions->sol_size_) + 1) <= storing_param->max_sol_size_)
        {
            return_code = store_solution_step_to_memory(
                system,
                simulation_status,
                solutions
            );
        }
        else
        {   
            return_code = ERROR_SOL_SIZE_EXCEED_MEMORY_ALLOC;
        }
    }
    else if (storing_param->storing_method_flag_ == STORING_METHOD_FLUSH)
    {
        return_code = flush_solution_step_to_csv_file(
            storing_param->flush_file_,
            system,
            simulation_status,
            solutions
        );
    }
    else if (storing_param->storing_method_flag_ == STORING_METHOD_DISABLED)
    {
        return_code = ERROR_STORE_SOLUTION_STEP_METHOD_DISABLED;
    }
    else
    {
        return_code = ERROR_STORE_SOLUTION_STEP_UNKNOWN_METHOD;
    }

    return return_code;
}

WIN32DLL_API int extend_sol_memory_buffer(
    Solutions *restrict solutions,
    StoringParam *restrict storing_param,
    const int objects_count
)
{
    int return_code;

    int64 buffer_size = storing_param->max_sol_size_ * 2;
    double *restrict temp_sol_state = NULL;
    double *restrict temp_sol_time = NULL;
    double *restrict temp_sol_dt = NULL;

    temp_sol_state = realloc(
        solutions->sol_state,
        buffer_size * objects_count * 6 * sizeof(double)
    );
    temp_sol_time = realloc(
        solutions->sol_time,
        buffer_size * sizeof(double)
    );
    temp_sol_dt = realloc(
        solutions->sol_dt,
        buffer_size * sizeof(double)
    );

    if (!temp_sol_state || !temp_sol_time || !temp_sol_dt)
    {
        return_code = ERROR_SOL_OUTPUT_EXTEND_MEMORY_REALLOC;
        goto error_memory;
    }

    storing_param->max_sol_size_ = buffer_size;
    solutions->sol_state = temp_sol_state;
    solutions->sol_time = temp_sol_time;
    solutions->sol_dt = temp_sol_dt;

    return SUCCESS;

error_memory:
    free(temp_sol_state);
    free(temp_sol_time);
    free(temp_sol_dt);
    return return_code;
}
