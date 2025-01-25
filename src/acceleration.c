#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "acceleration.h"
#include "acceleration_barnes_hut.h"
#include "acceleration_fast_multipole.h"

WIN32DLL_API int get_acceleration_method_flag(
    const char *restrict acceleration_method
)
{
    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        return ACCELERATION_METHOD_PAIRWISE;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        return ACCELERATION_METHOD_MASSLESS;
    }
    else if (strcmp(acceleration_method, "barnes_hut") == 0)
    {
        return ACCELERATION_METHOD_BARNES_HUT;
    }
    else if (strcmp(acceleration_method, "fast_multipole") == 0)
    {
        return ACCELERATION_METHOD_FAST_MULTIPOLE;
    }
    else
    {
        return -1;
    }
}

WIN32DLL_API void acceleration(
    int acceleration_method_flag,
    int objects_count,
    real *restrict x,
    real *restrict v,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length,
    real barnes_hut_theta
)
{
    (void) v;   // To suppress the warning of unused variable
    switch (acceleration_method_flag)
    {
        // Pairwise acceleration
        case ACCELERATION_METHOD_PAIRWISE:
            acceleration_pairwise(objects_count, x, a, m, G, softening_length);
            break;
        // Massless acceleration
        case ACCELERATION_METHOD_MASSLESS:
            acceleration_massless(objects_count, x, a, m, G, softening_length);
            break;
        // Barnes-Hut acceleration
        case ACCELERATION_METHOD_BARNES_HUT:
            acceleration_barnes_hut(objects_count, x, a, m, G, softening_length, barnes_hut_theta);
            break;
        case ACCELERATION_METHOD_FAST_MULTIPOLE:
            acceleration_fast_multipole(objects_count, x, a, m, G, softening_length);
            break;
        #ifdef USE_CUDA
            // CUDA Pairwise acceleration
            case ACCELERATION_METHOD_PAIRWISE_CUDA:
                acceleration_pairwise_cuda(objects_count, x, a, m, G, softening_length);
                break;
            // CUDA Pairwise acceleration with single precision
            case ACCELERATION_METHOD_PAIRWISE_FLOAT_CUDA:
                acceleration_pairwise_float_cuda(objects_count, x, a, m, G, softening_length);
                break;
        #endif
        default:
            fprintf(stderr, "Warning: Invalid acceleration method detected. Ignoring acceleration calculation.\n");
            break;
    }
}

WIN32DLL_API void acceleration_pairwise(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length
)
{   
    real R_norm;
    real temp_value;
    real temp_vec[3];
    real R[3];
    
    // Empty the input array
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    for(int i = 0; i < objects_count; i++)
    {
        for(int j = i + 1; j < objects_count; j++)
        {
            // Calculate \vec{R} and its norm
            R[0] = x[i * 3 + 0] - x[j * 3 + 0];
            R[1] = x[i * 3 + 1] - x[j * 3 + 1];
            R[2] = x[i * 3 + 2] - x[j * 3 + 2];
            R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2] + softening_length * softening_length);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
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
}

WIN32DLL_API void acceleration_massless(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length
)
{   
    real R_norm;
    real temp_value;
    real temp_vec[3];
    real R[3];
    
    // Empty the input array
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    int *restrict massive_indices = calloc(objects_count, sizeof(int));
    int *restrict massless_indices = calloc(objects_count, sizeof(int));
    int massive_objects_count = 0;
    int massless_objects_count = 0;
    for (int i = 0; i < objects_count; i++)
    {
        if (m[i] != 0.0)
        {
            massive_indices[massive_objects_count] = i;
            massive_objects_count++;
        }
        else
        {
            massless_indices[massless_objects_count] = i;
            massless_objects_count++;
        }
    }

    // Pairwise acceleration calculation for massive objects
    for (int i = 0; i < massive_objects_count; i++)
    {
        for (int j = i + 1; j < massive_objects_count; j++)
        {
            int idx_i = massive_indices[i];
            int idx_j = massive_indices[j];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2] + softening_length * softening_length);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[idx_i * 3 + 0] -= temp_vec[0] * m[j];
            a[idx_i * 3 + 1] -= temp_vec[1] * m[j];
            a[idx_i * 3 + 2] -= temp_vec[2] * m[j];
            a[idx_j * 3 + 0] += temp_vec[0] * m[i];
            a[idx_j * 3 + 1] += temp_vec[1] * m[i];
            a[idx_j * 3 + 2] += temp_vec[2] * m[i];
        }
    }

    // Acceleration calculation for massless objects
    for (int i = 0; i < massive_objects_count; i++)
    {
        for (int j = 0; j < massless_objects_count; j++)
        {
            int idx_i = massive_indices[i];
            int idx_j = massless_indices[j];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2] + softening_length * softening_length);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
            a[idx_j * 3 + 0] += temp_value * R[0] * m[i];
            a[idx_j * 3 + 1] += temp_value * R[1] * m[i];
            a[idx_j * 3 + 2] += temp_value * R[2] * m[i];
        }
    }

    free(massive_indices);
    free(massless_indices);
}
