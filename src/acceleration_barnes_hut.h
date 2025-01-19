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

#endif