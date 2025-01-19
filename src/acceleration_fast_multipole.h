#ifndef ACCELERATION_FAST_MULTIPOLE_H
#define ACCELERATION_FAST_MULTIPOLE_H

/**
 * acceleration_fast_multipole: Contains the fast multipole method 
 * for calculating gravitational acceleration
 */

#include "common.h"

void acceleration_fast_multipole(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length
);

#endif