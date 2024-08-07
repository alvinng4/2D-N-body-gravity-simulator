#ifndef COMMON_H
#define COMMON_H

/**
 * common: Contains commonly used functions for gravity 
 *         simulation, e.g. acceleration, initialize_system.
 *         Most definitions are also included in this file,
 *         including NPTS, int64, etc.
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#ifdef WIN32DLL_EXPORTS
    #define WIN32DLL_API __declspec(dllexport)
#else
    #define WIN32DLL_API 
#endif

#define BUFFER_SIZE 50000

typedef int64_t int64;
typedef double real;

typedef struct Solutions 
{
    double* sol_state;
    double* sol_time;
    double* sol_dt;
} Solutions;


/**
 * \brief Find the max absolute value in a 1D array
 * 
 * \param vec A 1D array
 * \param vec_length Length of the 1D array
 */
real abs_max_vec(const real *restrict vec, int vec_length);

/**
 * \brief Find the norm of a 1D array
 * 
 * \param vec A 1D array
 * \param vec_length Length of the 1D array
 */
real vec_norm(const real *restrict vec, int vec_length);

/**
 * \brief Compute the dot product of two 1D arrays
 * 
 * \param vec_1 A 1D array
 * \param vec_2 A 1D array
 * \param vec_length Length of the 1D arrays
 */
real vec_dot(
    const real *restrict vec_1,
    const real *restrict vec_2,
    int vec_length
);

/**
 * \brief Compute the cross product of two 1D arrays
 * 
 * \param vec_1 A 1D array
 * \param vec_2 A 1D array
 * \param result A 1D array to store the result
 */
void vec_cross(
    const real *restrict vec_1,
    const real *restrict vec_2,
    real *restrict result
);

/**
 * \brief Store the state of the system at a given time
 * 
 * \param file File pointer to the csv file
 * \param time time
 * \param dt Time step
 * \param x Array of position vectors of all objects
 * \param v Array of velocity vectors of all objects
 */
void write_to_csv_file(
    FILE *restrict file,
    double time,
    double dt,
    int objects_count,
    const double *restrict x,
    const double *restrict v,
    const double *restrict m,
    real G
);

/**
 * \brief Initialize default systems for gravity simulator.
 *        if the system name is recognized to be one of the 
 *        default system *x, *v and *m, *objects_count and 
 *        *G would be modified. Specifically, *x, *v and *m
 *        would be assigned a new block of memory.
 * 
 * \param system Name of the system to be initialized
 * \param x Pointer to pointer x, where x is the array of position vectors of all objects
 * \param v Pointer to pointer v, where v is the array of velocity vectors of all objects 
 * \param m Pointer to pointer m, where a is the array of acceleration vectors of all objects
 * \param objects_count Number of objects in the system
 * \param G Pointer to gravitational constant
 * 
 * \retval 0, if the system is successfully initialized
 * \retval 1, if the system is not recognized
 * \retval 2, if failed to allocate memory for x, v and m.
 */  
int initialize_system(
    const char *restrict system,
    real **x,
    real **v,
    real **m,
    int *restrict objects_count,
    real *restrict G
);

#endif