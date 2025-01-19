#ifndef ACCELERATION_BARNES_HUT_H
#define ACCELERATION_BARNES_HUT_H

/**
 * acceleration_barnes_hut: Contains the barnes-hut 
 * algorithm for calculating gravitational acceleration
 */

#include "common.h"

typedef struct BarnesHutTreeNode
{
    int index;
    real center_of_mass[3];
    real total_mass;
    real box_width;
    int internal_nodes_objects_count;
    struct BarnesHutTreeNode *children[8];
} BarnesHutTreeNode;

typedef struct BarnesHutTreeNodePool
{
    BarnesHutTreeNode *node_pool;
    int pool_size;
    struct BarnesHutTreeNodePool *next;
} BarnesHutTreeNodePool;    

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
    BarnesHutTreeNodePool *leaf_node_pool,
    BarnesHutTreeNodePool *internal_node_pool,
    int *restrict actual_interval_nodes_count,
    BarnesHutTreeNode *restrict root
);

// Shorten the tree if having less than 8 leaves
int _barnes_hut_shorten_tree(
    int *restrict actual_interval_nodes_count,
    BarnesHutTreeNode *restrict root
);

int _barnes_hut_compute_center_of_mass(
    int actual_interval_nodes_count,
    int *restrict max_depth,
    BarnesHutTreeNode *restrict root
);

int _barnes_hut_acceleration(
    real *restrict a,
    real G,
    real softening_length,
    real theta,
    int actual_interval_nodes_count,
    int max_depth,
    BarnesHutTreeNode *restrict root
);

void _barnes_hut_helper_acceleration_pair(
    BarnesHutTreeNode *restrict current_acc_leaf,
    BarnesHutTreeNode *restrict current_obj_leaf,
    real *restrict a,
    real G,
    real *restrict R,
    real R_norm
);

#endif