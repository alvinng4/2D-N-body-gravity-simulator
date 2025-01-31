/**
 * \file acceleration.h
 * \author Ching Yin Ng
 * \brief Contains simple methods in calculating gravitational
 *        acceleration, and acceleration-related functions.
 */

#ifndef ACCELERATION_H
#define ACCELERATION_H

#include "gravity_sim.h"

#define ACCELERATION_METHOD_PAIRWISE 0
#define ACCELERATION_METHOD_MASSLESS 1
#define ACCELERATION_METHOD_BARNES_HUT 2
#define ACCELERATION_METHOD_FAST_MULTIPOLE 3

/**
 * \brief Return acceleration method flag based on the input string
 * 
 * \param acceleration_method Name of the acceleration method
 * 
 * \retval SUCCESS If the acceleration method is recognized
 * \retval ERROR_UNKNOWN_ACCELERATION_METHOD If the acceleration method is not recognized
 */
int get_acceleration_method_flag(
    const char *restrict acceleration_method,
    uint *restrict acceleration_method_flag
);

/**
 * \brief Wrapper function for computing acceleration
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \retval SUCCESS If the computation is successful
 * \retval ERROR_UNKNOWN_ACCELERATION_CODE If the acceleration code is not recognized
 */
int acceleration(
    real *restrict a,
    const System *restrict system,
    AccelerationParam *restrict acceleration_param
);

#endif
