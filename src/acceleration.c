#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acceleration.h"
// #include "acceleration_barnes_hut.h"
// #include "acceleration_fast_multipole.h"
#include "error.h"
#include "gravity_sim.h"

/**
 * \brief Pairwise acceleration computation based on Newton's law of gravitational
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \retval SUCCESS If the computation is successful
 */
IN_FILE int acceleration_pairwise(
    real *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
);

/**
 * \brief Pairwise acceleration computation based on Newton's law of gravitational,
 *        ignoring the contribution of massless particles
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \retval SUCCESS If the computation is successful
 * \retval ERROR_ACCELERATION_MASSLESS_MEMORY_ALLOC If failed to allocate memory
 */
IN_FILE int acceleration_massless(
    real *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
);

WIN32DLL_API int get_acceleration_method_flag(
    const char *restrict acceleration_method,
    uint *restrict acceleration_method_flag
)
{
    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_PAIRWISE;
        return SUCCESS;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_MASSLESS;
        return SUCCESS;
    }
    else if (strcmp(acceleration_method, "barnes_hut") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_BARNES_HUT;
        return SUCCESS;
    }
    else if (strcmp(acceleration_method, "fast_multipole") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_FAST_MULTIPOLE;
        return SUCCESS;
    }
    else
    {
        return ERROR_UNKNOWN_ACCELERATION_METHOD;
    }
}

WIN32DLL_API int acceleration(
    real *restrict a,
    const System *restrict system,
    AccelerationParam *restrict acceleration_param
)
{
    switch (acceleration_param->acceleration_method_flag_)
    {
        case ACCELERATION_METHOD_PAIRWISE:
            return acceleration_pairwise(a, system, acceleration_param);
        case ACCELERATION_METHOD_MASSLESS:
            return acceleration_massless(a, system, acceleration_param);
        default:
            return ERROR_UNKNOWN_ACCELERATION_CODE;
    }
}

IN_FILE int acceleration_pairwise(
    real *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
)
{
    const int objects_count = system->objects_count;
    const real *x = system->x;
    const real *m = system->m;
    const real G = system->G;
    const real softening_length = acceleration_param->softening_length;

    /* Empty the input array */
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Compute the pairwise acceleration */
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = i + 1; j < objects_count; j++)
        {
            real temp_vec[3];
            real R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[i * 3 + 0] - x[j * 3 + 0];
            R[1] = x[i * 3 + 1] - x[j * 3 + 1];
            R[2] = x[i * 3 + 2] - x[j * 3 + 2];
            real R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            real temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[i * 3 + 0] -= temp_vec[0] * m[j];
            a[i * 3 + 1] -= temp_vec[1] * m[j];
            a[i * 3 + 2] -= temp_vec[2] * m[j];
            a[j * 3 + 0] += temp_vec[0] * m[i];
            a[j * 3 + 1] += temp_vec[1] * m[i];
            a[j * 3 + 2] += temp_vec[2] * m[i];
        }
    }

    return SUCCESS;
}

IN_FILE int acceleration_massless(
    real *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
)
{
    const int objects_count = system->objects_count;
    const real *x = system->x;
    const real *m = system->m;
    const real G = system->G;
    const real softening_length = acceleration_param->softening_length;

    /* Empty the input array */
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Find the numbers of massive and massless objects */
    int massive_objects_count = 0;
    int massless_objects_count = 0;
    for (int i = 0; i < objects_count; i++)
    {
        if (m[i] != 0.0)
        {
            massive_objects_count++;
        }
        else
        {
            massless_objects_count++;
        }
    }

    /* Find the indices of massive and massless objects */
    int *restrict massive_indices = malloc(massive_objects_count * sizeof(int));
    int *restrict massless_indices = malloc(massless_objects_count * sizeof(int));
    massive_objects_count = 0;
    massless_objects_count = 0;

    if (massive_indices == NULL || massless_indices == NULL)
    {
        goto malloc_error;
    }

    for (int i = 0; i < objects_count; i++)
    {
        if (m[i] != 0.0)
        {
            massless_indices[massless_objects_count] = i;
            massless_objects_count++;
        }
        else
        {
            massive_indices[massive_objects_count] = i;
            massive_objects_count++;
        }
    }

    /* Pairwise acceleration calculation for massive objects */
    for (int i = 0; i < massive_objects_count; i++)
    {
        for (int j = i + 1; j < massive_objects_count; j++)
        {
            int idx_i = massive_indices[i];
            int idx_j = massive_indices[j];
            real temp_vec[3];
            real R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            real R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            real temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[idx_i * 3 + 0] -= temp_vec[0] * m[idx_j];
            a[idx_i * 3 + 1] -= temp_vec[1] * m[idx_j];
            a[idx_i * 3 + 2] -= temp_vec[2] * m[idx_j];
            a[idx_j * 3 + 0] += temp_vec[0] * m[idx_i];
            a[idx_j * 3 + 1] += temp_vec[1] * m[idx_i];
            a[idx_j * 3 + 2] += temp_vec[2] * m[idx_i];
        }
    }

    /* Acceleration calculation for massless objects due to massive objects */
    for (int i = 0; i < massive_objects_count; i++)
    {
        for (int j = 0; j < massless_objects_count; j++)
        {
            int idx_i = massive_indices[i];
            int idx_j = massless_indices[j];
            real R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            real R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            real temp_value = G / (R_norm * R_norm * R_norm);
            a[idx_j * 3 + 0] += temp_value * R[0] * m[i];
            a[idx_j * 3 + 1] += temp_value * R[1] * m[i];
            a[idx_j * 3 + 2] += temp_value * R[2] * m[i];
        }
    }

    free(massive_indices);
    free(massless_indices);

    return SUCCESS;

malloc_error:
    free(massive_indices);
    free(massless_indices);
    return ERROR_ACCELERATION_MASSLESS_MEMORY_ALLOC;
}
