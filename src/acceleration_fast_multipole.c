#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "acceleration_fast_multipole.h"

WIN32DLL_API void acceleration_fast_multipole(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length
)
{       
    // Empty the input array
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    (void) x;
    (void) m;
    (void) G;
    (void) softening_length;
    fprintf(stderr, "Error: The fast multipole method is not implemented yet.\n");
}
