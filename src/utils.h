/**
 * \file utils.h
 * \brief Header file for utility functions
 * \author Ching Yin Ng
 */

#ifndef UTILS_H
#define UTILS_H

#include "gravity_sim.h"

/**
 * \brief Free memory allocated for a real array
 * 
 * \param ptr Pointer to the real array
 */
void free_memory_real(real *ptr);

/**
 * \brief Compute the energy of the system at a time step
 * 
 * \param system Pointer to the gravitational system
 * \param energy Pointer to the energy variable to be updated
 * 
 * \return SUCCESS if the energy is computed successfully
 */
int compute_energy_step(
    const System *restrict system,
    real *restrict energy
);

/**
 * \brief Compute the energy from solution state
 * 
 * \param objects_count Number of objects in the system
 * \param m Pointer to the mass array
 * \param G Gravitational constant
 * \param npts Number of time steps
 * \param count Pointer to the count variable
 * \param energy Pointer to the energy array to be updated
 * \param sol_state Pointer to the solution state array
 * \param is_exit Pointer to the exit flag
 */
void compute_energy_python(
    const int objects_count,
    const double *restrict m,
    const real G,
    const int npts,
    int *restrict count,
    real *restrict energy,
    const double (*restrict sol_state)[objects_count * 6],
    int *restrict is_exit
);

/**
 * \brief Compute the linear momentum from solution state
 * 
 * \param objects_count Number of objects in the system
 * \param m Pointer to the mass array
 * \param npts Number of time steps
 * \param count Pointer to the count variable
 * \param linear_momentum Pointer to the linear momentum array to be updated
 * \param sol_state Pointer to the solution state array
 * \param is_exit Pointer to the exit flag
 */
void compute_linear_momentum_python(
    const int objects_count,
    const double *restrict m,
    const int npts,
    int *restrict count,
    double *restrict linear_momentum,
    const double (*restrict sol_state)[objects_count * 6],
    int *restrict is_exit
);

/**
 * \brief Compute the angular momentum from solution state
 * 
 * \param objects_count Number of objects in the system
 * \param m Pointer to the mass array
 * \param npts Number of time steps
 * \param count Pointer to the count variable
 * \param angular_momentum Pointer to the angular momentum array to be updated
 * \param sol_state Pointer to the solution state array
 * \param is_exit Pointer to the exit flag
 */
void compute_angular_momentum_python(
    const int objects_count,
    const double *restrict m,
    const int npts,
    int *restrict count,
    double *restrict angular_momentum,
    const double (*restrict sol_state)[objects_count * 6],
    int *restrict is_exit
);

#endif
