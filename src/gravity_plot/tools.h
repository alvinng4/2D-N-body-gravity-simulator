#ifndef TOOLS_H
#define TOOLS_H

/**
 * tools: This file contain useful tools,
 *        such as compute_energy and free_memory_real.
 */

#include "common.h"

/**
 * \brief Frees allocated memory for pointers of type 'real'.
 * 
 *        This function should be called to free the allocated 
 *        memories after the python wrapper made a copy of the 
 *        arrays.
 * 
 * \param ptr Pointer to the memory block of type 'real' to be freed.
 * 
 * \return None
 */
void free_memory_real(real *ptr);

/**
 * \brief Compute the total energy in the newtonian system
 * 
 * \param objects_count Number of objects in the system
 * \param npts Length of the solution
 * \param count Pointer to the count of current progress
 * \param energy 1D array of size npts. Computed energy will 
 *               be stored into this array.
 * \param sol_state Array of state vectors of the solution
 * \param m Array of masses of the system
 * \param G Gravitational Constant
 * 
 * \return None
 */
void compute_energy(
    int objects_count, 
    int npts,
    int *restrict count, 
    double *restrict energy, 
    const double (*restrict sol_state)[objects_count * 6], 
    const double *restrict m, 
    real G
);

#endif