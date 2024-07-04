#include <math.h>
#include <stdlib.h>

#include "common.h"


WIN32DLL_API real abs_max_vec(const real *restrict vec, int vec_length)
{
    // Find the max absolute value in a 1D array
    real max = fabs(vec[0]);
    for (int i = 1; i < vec_length; i++)
    {
        max = fmax(max, fabs(vec[i]));
    }

    return max;
}

WIN32DLL_API real abs_max_vec_array(real (*restrict arr)[3], int objects_count)
{
    // Find the max absolute value in a 1D array
    real max = 0.0;
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            max = fmax(max, fabs(arr[i][j]));
        }
    }

    return max;
}

WIN32DLL_API real vec_norm(const real *restrict vec, int vec_length)
{   
    real sum = 0.0;
    if (vec_length == 3) 
    {
        sum = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    }
    else
    {
        for (int i = 0; i < vec_length; i++) sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}



WIN32DLL_API void acceleration(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict a)[3], 
    const real *restrict m, 
    real G
)
{   
    real R_norm, temp_value, temp_vec[3], R[3];

    // Empty the input array
    for (int i = 0; i < objects_count; i++)
    {
        a[i][0] = 0.0;
        a[i][1] = 0.0;
        a[i][2] = 0.0;
    }

    for(int i = 0; i < objects_count; i++)
    {
        for(int j = i + 1; j < objects_count; j++)
        {
            // Calculate \vec{R} and its norm
            R[0] = x[i][0] - x[j][0];
            R[1] = x[i][1] - x[j][1];
            R[2] = x[i][2] - x[j][2];
            R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[i][0] -= temp_vec[0] * m[j];
            a[i][1] -= temp_vec[1] * m[j];
            a[i][2] -= temp_vec[2] * m[j];
            a[j][0] += temp_vec[0] * m[i];
            a[j][1] += temp_vec[1] * m[i];
            a[j][2] += temp_vec[2] * m[i];
        }
    }
}
