#ifndef ACCELERATION_H
#define ACCELERATION_H

/**
 * acceleration: Contains methods in calculating acceleration
 */

#include "common.h"

#define ACCELERATION_METHOD_PAIRWISE 0
#define ACCELERATION_METHOD_MASSLESS 1
#define ACCELERATION_METHOD_BARNES_HUT 2
#define ACCELERATION_METHOD_FAST_MULTIPOLE 3
#ifdef USE_CUDA
    #define ACCELERATION_METHOD_PAIRWISE_CUDA 4
    #define ACCELERATION_METHOD_PAIRWISE_FLOAT_CUDA 5
#endif

/**
 * \brief Return acceleration method flag based on the input string
 * 
 * \param acceleration_method Acceleration method
 * 
 * \return Flag for acceleration method
 */
int get_acceleration_method_flag(
    const char *restrict acceleration_method
);

/**
 * \brief Wrapper function for acceleration computation
 * 
 * \param acceleration_method_flag Method to calculate acceleration (int flag)
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param v Array of velocity vectors of all objects
 * \param a Array of acceleration vectors to be modifed
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param barnes_hut_theta Theta parameter for Barnes-Hut algorithm
 */
void acceleration(
    int acceleration_method_flag,
    int objects_count,
    real *restrict x,
    real *restrict v,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length,
    real barnes_hut_theta
);

/**
 * \brief Pairwise computation of acceleration based on Newton's law of gravitational. 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param a Array of acceleration vectors to be modifed
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * 
 * \return None
 */
void acceleration_pairwise(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length
);

/**
 * \brief Compute acceleration based on Newton's law of gravitational,
 *        separating the calculation of massive and massless objects.
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param a Array of acceleration vectors to be modifed
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * 
 * \return None
 */
void acceleration_massless(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length
);

#endif
