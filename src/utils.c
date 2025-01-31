/**
 * \file utils.c
 * \brief Definitions of utility functions
 * \author Ching Yin Ng
 */
#include <stdlib.h>

#include "gravity_sim.h"
#include "math_functions.h"
#include "error.h"

WIN32DLL_API void free_memory_real(real *ptr)
{
    free(ptr);
}

WIN32DLL_API int compute_energy_step(
    const System *restrict system,
    real *restrict energy
)
{
    const real *restrict x = system->x;
    const real *restrict v = system->v;
    const real *restrict m = system->m;
    const real G = system->G;
    const int objects_count = system->objects_count;

    *energy = 0.0;

    for (int i = 0; i < objects_count; i++)
    {
        // KE
        real v_norm = vec_norm_3d(&v[i * 3]);
        *energy += (
            0.5 * m[i] 
            * v_norm * v_norm
        );

        // PE
        for (int j = i + 1; j < objects_count; j++)
        {
            real r_ij[3];
            for (int k = 0; k < 3; k++)
            {
                r_ij[k] = (x[i * 3 + k] - x[j * 3 + k]);
            }
            *energy -= (G * m[i] * m[j] / vec_norm_3d(r_ij));
        }
    }
    
    return SUCCESS;
}

WIN32DLL_API void compute_energy_python(
    const int objects_count,
    const double *restrict m,
    const real G,
    const int npts,
    int *restrict count,
    real *restrict energy,
    const double (*restrict sol_state)[objects_count * 6],
    int *restrict is_exit
)
{
    while (*count < npts)
    {   
        for (int i = 0; i < objects_count; i++)
        {
            // KE
            real v_norm = vec_norm_3d(&sol_state[*count][(objects_count + i) * 3]);
            energy[*count] += (
                0.5 * m[i] 
                * v_norm * v_norm
            );

            // PE
            for (int j = i + 1; j < objects_count; j++)
            {
                real r_ij[3];
                for (int k = 0; k < 3; k++)
                {
                    r_ij[k] = (
                        sol_state[*count][i * 3 + k]
                        - sol_state[*count][j * 3 + k]
                    );
                }
                energy[*count] -= (
                    G * m[i] * m[j]
                    / vec_norm_3d(r_ij)
                );
            }
        }
        *count += 1;

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            return;
        }
    }
}

WIN32DLL_API void compute_linear_momentum_python(
    const int objects_count,
    const double *restrict m,
    const int npts,
    int *restrict count,
    double *restrict linear_momentum,
    const double (*restrict sol_state)[objects_count * 6],
    int *restrict is_exit
)
{
    while (*count < npts)
    {
        real linear_momentum_vec_step[3] = {0.0};
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                linear_momentum_vec_step[j] += m[i] * (sol_state[*count][(objects_count + i) * 3 + j]);
            }
        }
        linear_momentum[*count] = vec_norm_3d(linear_momentum_vec_step);
        *count += 1;

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            return;
        }
    }
}

WIN32DLL_API void compute_angular_momentum_python(
    const int objects_count,
    const double *restrict m,
    const int npts,
    int *restrict count,
    double *restrict angular_momentum,
    const double (*restrict sol_state)[objects_count * 6],
    int *restrict is_exit
)
{
    while (*count < npts)
    {
        real angular_momentum_vec_step[3] = {0.0};
        // L = m * r x v
        for (int i = 0; i < objects_count; i++)
        {
            angular_momentum_vec_step[0] += m[i] * (
                sol_state[*count][i * 3 + 1] 
                * sol_state[*count][(objects_count + i) * 3 + 2]
                - sol_state[*count][i * 3 + 2] 
                * sol_state[*count][(objects_count + i) * 3 + 1]
            );
            angular_momentum_vec_step[1] += m[i] * (
                sol_state[*count][i * 3 + 2]
                * sol_state[*count][(objects_count + i) * 3 + 0]
                - sol_state[*count][i * 3 + 0]
                * sol_state[*count][(objects_count + i) * 3 + 2]
            );
            angular_momentum_vec_step[2] += m[i] * (
                sol_state[*count][i * 3]
                * sol_state[*count][(objects_count + i) * 3 + 1]
                - sol_state[*count][i * 3 + 1]
                * sol_state[*count][(objects_count + i) * 3]
            );
        }
        angular_momentum[*count] = vec_norm_3d(angular_momentum_vec_step);
        *count += 1;

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            return;
        }
    }
}
