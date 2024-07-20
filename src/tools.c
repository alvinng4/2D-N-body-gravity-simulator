#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "tools.h"

WIN32DLL_API void free_memory_real(real *ptr)
{
    free(ptr);
}

WIN32DLL_API void compute_energy(
    int objects_count, 
    int npts,
    int *restrict count, 
    double *restrict energy, 
    const double (*restrict sol_state)[objects_count * 6], 
    const double *restrict m, 
    real G,
    int *restrict is_exit
)
{
    real temp_vec[3];

    while (1)
    {   
        for (int i = 0; i < objects_count; i++)
        {
            // KE
            energy[*count] += (
                0.5 * m[i] 
                * pow(vec_norm(&sol_state[*count][(objects_count + i) * 3], 3), 2)
            );

            // PE
            for (int j = i + 1; j < objects_count; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    temp_vec[k] = (
                        sol_state[*count][i * 3 + k]
                        - sol_state[*count][j * 3 + k]
                    );
                }
                energy[*count] -= (
                    G * m[i] * m[j]
                    / vec_norm(temp_vec, 3)
                );
            }
        }
        *count += 1;

        if (*count >= npts)
        {
            break;
        }

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            return;
        }
    }
}


WIN32DLL_API void compute_linear_momentum(
    int objects_count, 
    int npts,
    int *restrict count,
    double *restrict linear_momentum,
    const double (*restrict sol_state)[objects_count * 6],
    const double *restrict m,
    int *restrict is_exit
)
{
    while (1)
    {   
        real temp_vec[3] = {0};
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                temp_vec[j] += m[i] * (sol_state[*count][(objects_count + i) * 3 + j]);
            }
        }
        linear_momentum[*count] = vec_norm(temp_vec, 3);
        *count += 1;

        if (*count >= npts)
        {
            break;
        }

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            return;
        }
    }
}


WIN32DLL_API void compute_angular_momentum(
    int objects_count, 
    int npts,
    int *restrict count,
    double *restrict angular_momentum,
    const double (*restrict sol_state)[objects_count * 6],
    const double *restrict m,
    int *restrict is_exit
)
{
    while (1)
    {   
        real temp_vec[3] = {0};
        for (int i = 0; i < objects_count; i++)
        {
            temp_vec[0] += m[i] * (
                sol_state[*count][i * 3 + 1]
                * sol_state[*count][(objects_count + i) * 3 + 2]
                - sol_state[*count][i * 3 + 2]
                * sol_state[*count][(objects_count + i) * 3 + 1]
            );
            temp_vec[1] += m[i] * (
                sol_state[*count][i * 3 + 2]
                * sol_state[*count][(objects_count + i) * 3 + 0]
                - sol_state[*count][i * 3 + 0]
                * sol_state[*count][(objects_count + i) * 3 + 2]
            );
            temp_vec[2] += m[i] * (
                sol_state[*count][i * 3]
                * sol_state[*count][(objects_count + i) * 3 + 1]
                - sol_state[*count][i * 3 + 1]
                * sol_state[*count][(objects_count + i) * 3]
            );
        }
        angular_momentum[*count] = vec_norm(temp_vec, 3);
        *count += 1;

        if (*count >= npts)
        {
            break;
        }

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            return;
        }
    }
}
