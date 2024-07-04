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
    real *restrict energy, 
    const real (*restrict sol_state)[objects_count * 6], 
    const real *restrict m, 
    real G
)
{
    // Round down current progress percentage as int
    int progress_percentage = (double) *count / npts * 100.0;

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

        // Exit to update progress bar
        if ((*count / npts * 100) > progress_percentage)
        {   
            break;
        }
    }
}
