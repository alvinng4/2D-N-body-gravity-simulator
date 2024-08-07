#ifndef ACCELERATION_H
#define ACCELERATION_H

/**
 * acceleration: Contains methods in calculating acceleration
 */

#include <stdbool.h>

#include "common.h"


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
    real G
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
    real G
);

/**
 * \brief Compute acceleration with Barnes-Hut algorithm
 * 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param a Array of acceleration vectors to be modifed
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param barnes_hut_theta Theta parameter for Barnes-Hut algorithm
 */
void acceleration_barnes_hut(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real barnes_hut_theta
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

#endif