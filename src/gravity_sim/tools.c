#include <math.h>

#include "common.h"

WIN32DLL_API real compute_energy(
    int objects_count, 
    const real (*restrict x)[3],
    const real (*restrict v)[3],
    const real *restrict m, 
    real G
)
{
    real temp_vec[3], energy = 0.0, norm;

    for (int i = 0; i < objects_count; i++)
    {   
        // KE
        norm = vec_norm(v[i], 3);
        if (norm != 0)
        {
            energy += (
                0.5 * m[i] 
                * pow(norm, 2)
            );
        }
        else
        {
            return NAN;
        }

        // PE
        for (int j = i + 1; j < objects_count; j++)
        {
            temp_vec[0] = (
                x[i][0] 
                - x[j][0]
            );
            temp_vec[1] = (
                x[i][1] 
                - x[j][1]
            );
            temp_vec[2] = (
                x[i][2] 
                - x[j][2]
            );

            norm = vec_norm(temp_vec, 3);
            if (norm != 0)
            {
                energy -= (
                    G * m[i] * m[j]
                    / norm
                );
            }
            else
            {
                return NAN;
            }
        }
    }

    return energy;
    
}