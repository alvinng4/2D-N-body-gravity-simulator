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
#define BARNES_HUT_THETA 0.5

typedef int64_t int64;
typedef double real;

typedef struct Solutions 
{
    double* sol_state;
    double* sol_time;
    double* sol_dt;
} Solutions;

typedef struct BarnesHutTreeNode
{
    bool is_leaf;
    int index;
    real total_mass;
    real center_of_mass[3];
    real box_width;
    struct BarnesHutTreeNode *children[8];
} BarnesHutTreeNode;

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
 * \brief Pairwise computation of acceleration based on Newton's law of gravitational. 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param a Array of acceleration vectors of all objects.
 *                   Values inside this array will be replaced
 *                   by the computed acceleration.
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
    real G
);

/**
 * \brief Compute acceleration based on Newton's law of gravitational,
 *        separating the calculation of massive and massless objects.
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param a Array of acceleration vectors of all objects.
 *                   Values inside this array will be replaced
 *                   by the computed acceleration.
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
    real G
);

void acceleration_barnes_hut(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G
);

int _barnes_hut_check_quadrant(
    real x,
    real y,
    real z,
    real center_x,
    real center_y,
    real center_z
);

int _barnes_hut_construct_octree(
    int objects_count,
    const real *restrict x,
    const real *restrict m,
    real width,
    BarnesHutTreeNode *root
);

int _barnes_hut_compute_center_of_mass(BarnesHutTreeNode *root);

int _barnes_hut_acceleration(
    real theta,
    int objects_count,
    real *restrict a,
    real G,
    BarnesHutTreeNode *root
);

void _barnes_hut_helper_acceleration_pair(
    BarnesHutTreeNode *current_acc_leaf,
    BarnesHutTreeNode *current_obj_leaf,
    real *restrict a,
    real G,
    real *restrict R,
    real R_norm
);

WIN32DLL_API int _barnes_hut_free_octree(BarnesHutTreeNode *restrict root);

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