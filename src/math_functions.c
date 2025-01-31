#include <math.h>

#include "gravity_sim.h"


WIN32DLL_API real vec_norm_3d(const real *restrict vec)
{
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

WIN32DLL_API real vec_norm(const real *restrict vec, int vec_length)
{   
    real sum = 0.0;
    for (int i = 0; i < vec_length; i++) 
    {
        sum += vec[i] * vec[i];
    }

    return sqrt(sum);
}

WIN32DLL_API real abs_max_vec(const real *restrict vec, int vec_length)
{
    real max = fabs(vec[0]);
    for (int i = 1; i < vec_length; i++)
    {
        max = fmax(max, fabs(vec[i]));
    }

    return max;
}

WIN32DLL_API real vec_dot_3d(
    const real *restrict vec_1,
    const real *restrict vec_2
)
{
    return vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2];
}

WIN32DLL_API real vec_dot(
    const real *restrict vec_1,
    const real *restrict vec_2,
    const int vec_length
)
{
    real sum = 0.0;
    for (int i = 0; i < vec_length; i++)
    {
        sum += vec_1[i] * vec_2[i];
    }

    return sum;
}