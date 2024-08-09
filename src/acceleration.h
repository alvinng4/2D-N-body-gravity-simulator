#ifndef ACCELERATION_H
#define ACCELERATION_H

/**
 * acceleration: Contains methods in calculating acceleration
 */

#include <stdbool.h>

#include "common.h"

#define ACCELERATION_METHOD_PAIRWISE 0
#define ACCELERATION_METHOD_MASSLESS 1
#define ACCELERATION_METHOD_BARNES_HUT 2

typedef struct BarnesHutTreeNode
{
    int index;
    real total_mass;
    real center_of_mass[3];
    real box_width;
    struct BarnesHutTreeNode *children[8];
} BarnesHutTreeNode;

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

/**
 * \brief Compute acceleration with Barnes-Hut algorithm
 * 
 * \param objects_count Number of objects
 * \param x Array of position vectors
 * \param a Array of acceleration vectors to be modifed
 * \param m Array of masses
 * \param G Gravitational constant
 * \param barnes_hut_theta Theta parameter for Barnes-Hut algorithm
 */
void acceleration_barnes_hut(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length,
    real barnes_hut_theta
);

/**
 * \brief Calculate the bounding box of the system
 * 
 * \param objects_count Number of objects
 * \param x Array of position vectors
 */
void _calculate_bounding_box(
    int objects_count,
    const real *restrict x,
    real *restrict center,
    real *restrict width
);

/**
 * \brief Check the region relative to the center
 *  
 * \param x x-coordinate
 * \param y y-coordinate
 * \param z z-coordinate
 * \param center_x x-coordinate of the center
 * \param center_y y-coordinate of the center
 * \param center_z z-coordinate of the center
 * 
 * \return Region index
 */
int _barnes_hut_check_region(
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
    BarnesHutTreeNode *restrict root
);

int _barnes_hut_compute_center_of_mass(BarnesHutTreeNode *root);

int _barnes_hut_acceleration(
    real theta,
    int objects_count,
    real *restrict a,
    real G,
    real softening_length,
    BarnesHutTreeNode *restrict root
);

void _barnes_hut_helper_acceleration_pair(
    BarnesHutTreeNode *restrict current_acc_leaf,
    BarnesHutTreeNode *restrict current_obj_leaf,
    real *restrict a,
    real G,
    real softening_length,
    real *restrict R,
    real R_norm
);

WIN32DLL_API int _barnes_hut_free_octree(BarnesHutTreeNode *restrict root);

#endif