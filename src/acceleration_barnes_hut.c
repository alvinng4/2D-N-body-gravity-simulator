#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "gravity_sim.h"
#include "math_functions.h"

#define MAX_NUM_PARTICLES_PER_LEAF 8

typedef struct BarnesHutTreeLeaf
{
    int objects_count;
    int indices[MAX_NUM_PARTICLES_PER_LEAF];
} BarnesHutTreeLeaf;

typedef struct BarnesHutTreeLeafPool
{
    BarnesHutTreeLeaf *leaves;
    int pool_size;
    struct BarnesHutTreeLeafPool *next;
} BarnesHutTreeLeafPool;

typedef struct BarnesHutTreeNode
{
    real center_of_mass[3];
    real total_mass;
    real box_width;
    struct BarnesHutTreeNode *children[8];
    BarnesHutTreeLeaf *leaves[8];
} BarnesHutTreeNode;

typedef struct BarnesHutTreeNodePool
{
    BarnesHutTreeNode *nodes;
    int pool_size;
    struct BarnesHutTreeNodePool *next;
} BarnesHutTreeNodePool; 

typedef struct BarnesHutCOMStack
{
    BarnesHutTreeNode *node;
    struct BarnesHutCOMStack *last;
    real sum_of_mass_times_distance[3];
    int processed_region;
} BarnesHutCOMStack;

typedef struct BarnesHutAccStack
{
    BarnesHutTreeNode *node;
    struct BarnesHutAccStack *last;
    int processed_region;
} BarnesHutAccStack;

// For debug
// Helper function to print indentation based on the tree level.
static void print_indent(int level) {
    for (int i = 0; i < level; i++) {
        printf("    ");
    }
}

// Recursively visualize the octree.
void visualize_octree(BarnesHutTreeNode *node, int level, real *x, real *m) {
    if (node == NULL)
        return;

    print_indent(level);
    printf("Node: center_of_mass=[%f, %f, %f], total_mass=%g, box_width=%f\n",
           node->center_of_mass[0], node->center_of_mass[1], node->center_of_mass[2],
           node->total_mass, node->box_width);

    // Iterate through each of the 8 octants.
    for (int i = 0; i < 8; i++) {
        if (node->children[i] != NULL) {
            print_indent(level);
            printf("Child[%d]:\n", i);
            visualize_octree(node->children[i], level + 1, x, m);
        }
        if (node->leaves[i] != NULL) {
            print_indent(level);
            printf("Leaf[%d]: objects_count = %d\n", i, node->leaves[i]->objects_count);
            for (int j = 0; j < node->leaves[i]->objects_count; j++)
            {
                print_indent(level + 1);
                printf(
                    "idx = %d, x = [%f, %f, %f], m = %g\n",
                       node->leaves[i]->indices[j],
                       x[node->leaves[i]->indices[j] * 3 + 0],
                       x[node->leaves[i]->indices[j] * 3 + 1],
                       x[node->leaves[i]->indices[j] * 3 + 2],
                       m[node->leaves[i]->indices[j]]
                    );
            }
        }
    }
}

/**
 * \brief Calculate the bounding box of the system
 * 
 * \param objects_count Number of objects
 * \param x Array of position vectors
 * \param center 3D vector of the center of the bounding box
 * \param width Width of the bounding box
 */
IN_FILE void _calculate_bounding_box(
    const int objects_count,
    const real *restrict x,
    real *restrict center,
    real *restrict width
)
{
    /* Find the width of the bounding box */
    real min_x = x[0];
    real max_x = x[0];
    real min_y = x[1];
    real max_y = x[1];
    real min_z = x[2];
    real max_z = x[2];

    for (int i = 1; i < objects_count; i++)
    {
        min_x = fmin(min_x, x[i * 3 + 0]);
        max_x = fmax(max_x, x[i * 3 + 0]);
        min_y = fmin(min_y, x[i * 3 + 1]);
        max_y = fmax(max_y, x[i * 3 + 1]);
        min_z = fmin(min_z, x[i * 3 + 2]);
        max_z = fmax(max_z, x[i * 3 + 2]);
    }

    center[0] = (max_x + min_x) / 2.0;
    center[1] = (max_y + min_y) / 2.0;
    center[2] = (max_z + min_z) / 2.0;

    real width_x = max_x - min_x;
    real width_y = max_y - min_y;
    real width_z = max_z - min_z;
    *width = fmax(fmax(width_x, width_y), width_z);
}

/**
 * \brief Initialize Barnes-Hut tree node
 * 
 * \param node Pointer to the node
 * \param center_of_mass 3d vector of the center of mass
 * \param total_mass Total mass
 * \param box_width Width of the box
 */
IN_FILE void _initialize_node(
    BarnesHutTreeNode *restrict node,
    const real center_of_mass[3],
    const real total_mass,
    const real box_width
)
{
    node->center_of_mass[0] = center_of_mass[0];
    node->center_of_mass[1] = center_of_mass[1];
    node->center_of_mass[2] = center_of_mass[2];
    node->total_mass = total_mass;
    node->box_width = box_width;
    for (int i = 0; i < 8; i++)
    {
        node->children[i] = NULL;
        node->leaves[i] = NULL;
    }
}

/**
 * \brief Check the region (i.e. octant)
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
IN_FILE int _check_region(
    const real x,
    const real y,
    const real z,
    const real center_x,
    const real center_y,
    const real center_z
)
{
    int region = 0;

    if (x > center_x)
    {
        region |= 1;
    }
    if (y > center_y)
    {
        region |= 2;
    }
    if (z > center_z) 
    {
        region |= 4;
    }

    return region;
}

/**
 * \brief Get leaf node from leaf pool
 * 
 * \param leaf Pointer of Pointer to the leaf node
 * \param leaf_pool_ptr Pointer to pool of leaf nodes
 * \param leaf_pool_expand_count Pointer to the number of leaf nodes after expanding the pool
 */
IN_FILE int _get_leaf(
    BarnesHutTreeLeaf **restrict leaf,
    BarnesHutTreeLeafPool **leaf_pool_ptr,
    int *restrict leaf_pool_expand_count
)
{
    int return_code;

    const int pool_size = (*leaf_pool_ptr)->pool_size;
    if (*leaf_pool_expand_count >= pool_size)
    {
        (*leaf_pool_ptr)->next = malloc(sizeof(BarnesHutTreeLeafPool));
        if (!((*leaf_pool_ptr)->next))
        {
            return_code = ERROR_BARNES_HUT_GET_LEAF_POOL_PTR_MEMORY_ALLOC;
            goto err_leaf_pool_ptr_memory;
        }
        *leaf_pool_ptr = (*leaf_pool_ptr)->next;
        (*leaf_pool_ptr)->next = NULL;
        (*leaf_pool_ptr)->pool_size = pool_size * 2;
        (*leaf_pool_ptr)->leaves = malloc((*leaf_pool_ptr)->pool_size * sizeof(BarnesHutTreeLeaf));
        if (!((*leaf_pool_ptr)->leaves))
        {
            return_code = ERROR_BARNES_HUT_GET_LEAF_POOL_MEMORY_ALLOC;
            goto err_leaf_pool_memory;
        }

        *leaf_pool_expand_count = 0;
    }

    *leaf = &((*leaf_pool_ptr)->leaves[*leaf_pool_expand_count]);
    (*leaf)->objects_count = 0;
    *leaf_pool_expand_count += 1;

    return SUCCESS;

// The memory will be freed in the main 
// function, even in the case of error
err_leaf_pool_memory:
err_leaf_pool_ptr_memory:
    return return_code;
}

/**
 * \brief Divide the tree when the leaf is full
 * 
 * \param node Pointer to the node
 * \param region Region index of the full leaf
 * \param node_pool_ptr Pointer to pool of nodes
 * \param node_pool_expand_count Pointer to the number of nodes after expanding the pool
 * \param leaf_pool_ptr Pointer to pool of leaves
 * \param leaf_pool_expand_count Pointer to the number of leaves after expanding the pool
 * \param x Array of position vectors
 * \param m Array of masses
 */
IN_FILE int _divide_tree(
    BarnesHutTreeNode *restrict const node,
    const int region,
    BarnesHutTreeNodePool **node_pool_ptr,
    int *restrict node_pool_expand_count,
    BarnesHutTreeLeafPool **leaf_pool_ptr,
    int *restrict leaf_pool_expand_count,
    const real *restrict x,
    const real *restrict m
)
{
    int return_code;

    int pool_size = (*node_pool_ptr)->pool_size;
    if (*node_pool_expand_count >= pool_size)
    {
        (*node_pool_ptr)->next = malloc(sizeof(BarnesHutTreeNodePool));
        if (!((*node_pool_ptr)->next))
        {
            return_code = ERROR_BARNES_HUT_DIVIDE_LEAF_NODE_POOL_PTR_MEMORY_ALLOC;
            goto err_node_pool_ptr_memory;
        }
        (*node_pool_ptr) = (*node_pool_ptr)->next;
        (*node_pool_ptr)->next = NULL;
        (*node_pool_ptr)->pool_size = pool_size * 2;
        (*node_pool_ptr)->nodes = malloc((*node_pool_ptr)->pool_size * sizeof(BarnesHutTreeNode));
        if ((*node_pool_ptr)->nodes == NULL)
        {
            return_code = ERROR_BARNES_HUT_DIVIDE_LEAF_NODE_POOL_MEMORY_ALLOC;
            goto err_node_pool_memory;
        }

        *node_pool_expand_count = 0;
    }

    BarnesHutTreeNode *new_node = &((*node_pool_ptr)->nodes[*node_pool_expand_count]);
    *node_pool_expand_count += 1;
    _initialize_node(new_node, node->center_of_mass, 0.0, node->box_width / 2.0);
    if (region & 1)
    {
        new_node->center_of_mass[0] += new_node->box_width / 2.0;
    }
    else
    {
        new_node->center_of_mass[0] -= new_node->box_width / 2.0;
    }

    if (region & 2)
    {
        new_node->center_of_mass[1] += new_node->box_width / 2.0;
    }
    else
    {
        new_node->center_of_mass[1] -= new_node->box_width / 2.0;
    }

    if (region & 4)
    {
        new_node->center_of_mass[2] += new_node->box_width / 2.0;
    }
    else
    {
        new_node->center_of_mass[2] -= new_node->box_width / 2.0;
    }


    BarnesHutTreeLeaf *old_leaf = node->leaves[region];
    const int num_objects_old_leaf = old_leaf->objects_count;
    int old_leaf_indices[MAX_NUM_PARTICLES_PER_LEAF];
    memcpy(old_leaf_indices, old_leaf->indices, num_objects_old_leaf * sizeof(int));

    // Reset the old leaf
    old_leaf->objects_count = 0;
    node->children[region] = new_node;
    node->leaves[region] = NULL;

    // Distribute the objects for the new node
    for (int i = 0; i < num_objects_old_leaf; i++)
    {
        const int idx_i = old_leaf_indices[i];
        const int new_region = _check_region(
            x[idx_i * 3 + 0],
            x[idx_i * 3 + 1],
            x[idx_i * 3 + 2],
            new_node->center_of_mass[0],
            new_node->center_of_mass[1],
            new_node->center_of_mass[2]
        );

        // Check if the leaf is not created
        BarnesHutTreeLeaf *new_leaf = new_node->leaves[new_region];
        if (!new_leaf)
        {
            if (i == 0)
            {
                new_leaf = old_leaf;
            }
            else
            {
                return_code = _get_leaf(
                    &new_leaf,
                    leaf_pool_ptr,
                    leaf_pool_expand_count
                );
                if (return_code != SUCCESS)
                {
                    goto err_get_leaf;
                }
            }
        }
        
        new_leaf->indices[new_leaf->objects_count] = idx_i;
        new_leaf->objects_count += 1;

        new_node->leaves[new_region] = new_leaf;
        new_node->total_mass += m[idx_i];
    }

    return SUCCESS;


// The memory will be freed in the main
// function, even in the case of error
err_get_leaf:
err_node_pool_ptr_memory:
err_node_pool_memory:
    return return_code;
}

/**
 * \brief Construct the octree
 * 
 * \param objects_count Number of objects
 * \param x Array of position vectors
 * \param m Array of masses
 * \param leaf_pool Pool of leaves
 * \param node_pool Pool of nodes
 * \param max_depth Pointer to the maximum depth of the tree
 * \param root Root node of the octree
 */
IN_FILE int _construct_octree(
    const int objects_count,
    const real *restrict x,
    const real *restrict m,
    BarnesHutTreeLeafPool *leaf_pool,
    BarnesHutTreeNodePool *node_pool,
    int *restrict max_depth,
    BarnesHutTreeNode *restrict root
)
{
    int return_code;

    // To keep track of the number of nodes / leaves
    // after expanding the pool
    int leaf_pool_expand_count = 0;
    int node_pool_expand_count = 0;

    // Construct octree
    for (int i = 0; i < objects_count; i++)
    {
        BarnesHutTreeNode *node = root;
        int depth = 0;

        while (true)
        {
            depth++;

            int region = _check_region(
                x[i * 3 + 0],
                x[i * 3 + 1],
                x[i * 3 + 2],
                node->center_of_mass[0],
                node->center_of_mass[1],
                node->center_of_mass[2]
            );

            BarnesHutTreeNode *child = node->children[region];
            BarnesHutTreeLeaf *leaf = node->leaves[region];

            node->total_mass += m[i];

            // Both child node and leaf do not exist
            if ((!child) && (!leaf))
            {
                // Create a new leaf
                return_code = _get_leaf(
                    &leaf,
                    &leaf_pool,
                    &leaf_pool_expand_count
                );
                if (return_code != SUCCESS)
                {
                    goto err_get_leaf;
                }
                leaf->indices[0] = i;
                leaf->objects_count = 1;
                node->leaves[region] = leaf;
                break;
            }

            // Leaf exists
            else if (!child)
            {
                // Check if the leaf is full
                if (leaf->objects_count >= MAX_NUM_PARTICLES_PER_LEAF)
                {
                    // Split the branch
                    return_code = _divide_tree(
                        node,
                        region,
                        &node_pool,
                        &node_pool_expand_count,
                        &leaf_pool,
                        &leaf_pool_expand_count,
                        x,
                        m
                    );
                    if (return_code != SUCCESS)
                    {
                        goto err_divide_leaf;
                    }
                    node = node->children[region];
                }
                else
                {
                    leaf->indices[leaf->objects_count] = i;
                    leaf->objects_count += 1;

                    break;
                }
            }

            // Child node already exists
            else
            {
                node = node->children[region];
            }
        }

        if (depth > *max_depth)
        {
            *max_depth = depth;
        }
    }

    return SUCCESS;

err_divide_leaf:
err_get_leaf:
    return return_code;
}

/**
 * \brief Update the sum of mass times distance vector
 * 
 * \param sum_of_mass_times_distance_vector Pointer to the sum of mass times distance vector
 * \param x Array of position vectors
 * \param m Array of masses
 * \param leaf Leaf node
 */
IN_FILE void _com_update_sum_of_mass_times_distance(
    real *restrict sum_of_mass_times_distance_vector,
    const real *restrict x,
    const real *restrict m,
    const BarnesHutTreeLeaf *restrict leaf
)
{
    for (int i = 0; i < (leaf->objects_count); i++)
    {
        const int idx_i = leaf->indices[i];
        real m_i = m[idx_i];
        sum_of_mass_times_distance_vector[0] += m_i * x[idx_i * 3 + 0];
        sum_of_mass_times_distance_vector[1] += m_i * x[idx_i * 3 + 1];
        sum_of_mass_times_distance_vector[2] += m_i * x[idx_i * 3 + 2];
    }
}

/**
 * \brief Compute the center of mass
 * 
 * \param x Array of position vectors
 * \param m Array of masses
 * \param max_depth Maximum depth of the tree
 * \param root Root node of the octree
 */
IN_FILE int _compute_center_of_mass(
    const real *restrict x,
    const real *restrict m,
    const int max_depth,
    BarnesHutTreeNode *restrict root
)
{
    int return_code;

    BarnesHutTreeNode *node = root;

    /* Create a stack pool */
    const int stack_pool_size = max_depth;
    BarnesHutCOMStack *stack_pool = malloc(stack_pool_size * sizeof(BarnesHutCOMStack));
    if (!stack_pool)
    {
        return_code = ERROR_BARNES_HUT_COMPUTE_COM_STACK_MEMORY_ALLOC;
        goto err_memory;
    }

    int stack_count = 1;
    BarnesHutCOMStack *stack = &(stack_pool[0]);
    stack->node = root;
    stack->last = NULL;
    stack->sum_of_mass_times_distance[0] = 0.0;
    stack->sum_of_mass_times_distance[1] = 0.0;
    stack->sum_of_mass_times_distance[2] = 0.0;
    stack->processed_region = -1;

    while (true)
    {
        for (int i = (stack->processed_region + 1); i < 8; i++)
        {
            BarnesHutTreeNode *child_i = node->children[i];
            BarnesHutTreeLeaf *leaf_i = node->leaves[i];

            // Both child node and leaf do not exist
            if ((!child_i) && (!leaf_i))
            {
                stack->processed_region = i;
            }

            // Leaf exists
            else if (!child_i)
            {
                _com_update_sum_of_mass_times_distance(
                    stack->sum_of_mass_times_distance,
                    x,
                    m,
                    leaf_i
                );
                stack->processed_region = i;
            }

            // Node exists
            else
            {
                // Create a new stack item
                if (stack_count >= stack_pool_size)
                {
                    return_code = ERROR_BARNES_HUT_COMPUTE_COM_STACK_FULL;
                    goto err_stack_full;
                }
                BarnesHutCOMStack *new_item = &stack_pool[stack_count];
                stack_count++;

                new_item->node = child_i;
                new_item->last = stack;
                new_item->sum_of_mass_times_distance[0] = 0.0;
                new_item->sum_of_mass_times_distance[1] = 0.0;
                new_item->sum_of_mass_times_distance[2] = 0.0;
                new_item->processed_region = -1;

                stack = new_item;
                node = child_i;

                break;
            }
        }

        if (stack->processed_region >= 7)
        {
            BarnesHutCOMStack *parent_stack = stack->last;

            // Break if the stack is empty
            if (!parent_stack)
            {
                break;
            }

            real total_mass = node->total_mass;
            node->center_of_mass[0] = stack->sum_of_mass_times_distance[0] / total_mass;
            node->center_of_mass[1] = stack->sum_of_mass_times_distance[1] / total_mass;
            node->center_of_mass[2] = stack->sum_of_mass_times_distance[2] / total_mass;

            parent_stack->sum_of_mass_times_distance[0] += stack->sum_of_mass_times_distance[0];
            parent_stack->sum_of_mass_times_distance[1] += stack->sum_of_mass_times_distance[1];
            parent_stack->sum_of_mass_times_distance[2] += stack->sum_of_mass_times_distance[2];

            parent_stack->processed_region += 1;
            stack = parent_stack;
            node = parent_stack->node;
            stack_count--;
        }
    }

    /* Free the stack pool */
    free(stack_pool);

    return SUCCESS;

err_stack_full:
err_memory:
    free(stack_pool);
    return return_code;
}

/**
 * \brief Check if the tree walk node includes the given leaf
 * 
 * \param branch_record_leaf Branch record of the given leaf
 * \param branch_record_leaf_count Size of branch record
 * \param branch_record_tree_walk Branch record of the tree walk node / leaf
 * \param branch_record_tree_walk_count Size of branch record
 * 
 * \return True if the tree walk node includes the given leaf
 * \return False otherwise
 */
IN_FILE bool _compare_branch_record(
    const int *restrict branch_record_leaf,
    const int branch_record_leaf_count,
    const int *restrict branch_record_tree_walk,
    const int branch_record_tree_walk_count
)
{
    if (branch_record_tree_walk_count > branch_record_leaf_count)
    {
        return false;
    }

    for (int i = 0; i < branch_record_tree_walk_count; i++)
    {
        if (branch_record_tree_walk[i] != branch_record_leaf[i])
        {
            return false;
        }
    }

    return true;
}

/**
 * \brief Compute the acceleration between all particles within the leaf
 * 
 * \param a Array of acceleration vectors to be modified
 * \param x Array of position vectors
 * \param m Array of masses
 * \param G Gravitational constant
 * \param softening_length Softening length
 * \param leaf Leaf node
 */
IN_FILE void _compute_acc_single_leaf(
    real *restrict a,
    const real *restrict x,
    const real *restrict m,
    const real G,
    const real softening_length,
    BarnesHutTreeLeaf *restrict leaf
)
{
    const int objects_count = leaf->objects_count;
    const int *restrict indices = leaf->indices;

    /* Compute the pairwise acceleration */
    for (int i = 0; i < objects_count; i++)
    {
        const int idx_i = indices[i];
        const real m_i = m[idx_i];
        for (int j = i + 1; j < objects_count; j++)
        {
            const int idx_j = indices[j];
            const real m_j = m[idx_j];

            real temp_vec[3];
            real R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            const real R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            const real temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[idx_i * 3 + 0] -= temp_vec[0] * m_j;
            a[idx_i * 3 + 1] -= temp_vec[1] * m_j;
            a[idx_i * 3 + 2] -= temp_vec[2] * m_j;
            a[idx_j * 3 + 0] += temp_vec[0] * m_i;
            a[idx_j * 3 + 1] += temp_vec[1] * m_i;
            a[idx_j * 3 + 2] += temp_vec[2] * m_i;
        }
    }
}

/**
 * \brief Compute the acceleration of particles of leaf_i due to leaf_j
 * 
 * \param a Array of acceleration vectors to be modified
 * \param x Array of position vectors
 * \param m Array of masses
 * \param G Gravitational constant
 * \param softening_length Softening length
 * \param indices_i Array of indices of particles in leaf_i
 * \param objects_count_i Number of particles in leaf_i
 * \param leaf Leaf node
 */
IN_FILE void _compute_acc_leaf_to_leaf(
    real *restrict a,
    const real *restrict x,
    const real *restrict m,
    const real G,
    const real softening_length,
    const int *restrict indices_i,
    const int objects_count_i,
    BarnesHutTreeLeaf *restrict leaf_j
)
{
    const int objects_count_j = leaf_j->objects_count;
    const int *restrict indices_j = leaf_j->indices;

    /* Compute the pairwise acceleration */
    for (int i = 0; i < objects_count_i; i++)
    {
        const int idx_i = indices_i[i];
        for (int j = 0; j < objects_count_j; j++)
        {
            const int idx_j = indices_j[j];
            const real m_j = m[idx_j];

            real temp_vec[3];
            real R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            const real R_norm = sqrt(
                R[0] * R[0] +
                R[1] * R[1] +
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            const real temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[idx_i * 3 + 0] -= temp_vec[0] * m_j;
            a[idx_i * 3 + 1] -= temp_vec[1] * m_j;
            a[idx_i * 3 + 2] -= temp_vec[2] * m_j;
        }
    }
}

// /**
//  * \brief Compute the acceleration of particles of leaf due to node
//  * 
//  * \param a Array of acceleration vectors to be modified
//  * \param x Array of position vectors
//  * \param G Gravitational constant
//  * \param softening_length Softening length
//  * \param leaf Leaf node
//  * \param center_of_mass Center of mass of the node
//  * \param m_j Total mass of the node
//  */
// IN_FILE void _compute_acc_leaf_to_node(
//     real *restrict a,
//     const real *restrict x,
//     const real G,
//     const real softening_length,
//     BarnesHutTreeLeaf *restrict leaf,
//     const real *restrict center_of_mass,
//     const real m_j
// )
// {
//     const int objects_count = leaf->objects_count;
//     const int *restrict indices = leaf->indices;

//     for (int i = 0; i < objects_count; i++)
//     {
//         real R[3];
//         real temp_vec[3];

//         const int idx_i = indices[i];

//         // Calculate \vec{R} and its norm
//         R[0] = x[idx_i * 3 + 0] - center_of_mass[0];
//         R[1] = x[idx_i * 3 + 1] - center_of_mass[1];
//         R[2] = x[idx_i * 3 + 2] - center_of_mass[2];
//         const real R_norm = sqrt(
//             R[0] * R[0] + 
//             R[1] * R[1] + 
//             R[2] * R[2] +
//             softening_length * softening_length
//         );

//         // Calculate the acceleration
//         const real temp_value = G / (R_norm * R_norm * R_norm);
//         temp_vec[0] = temp_value * R[0];
//         temp_vec[1] = temp_value * R[1];
//         temp_vec[2] = temp_value * R[2];
//         a[idx_i * 3 + 0] -= temp_vec[0] * m_j;
//         a[idx_i * 3 + 1] -= temp_vec[1] * m_j;
//         a[idx_i * 3 + 2] -= temp_vec[2] * m_j;
//     }
// }

IN_FILE void _pop_particles(
    int *restrict particle_indices,
    int *restrict num_particles,
    int *restrict pop_indices,
    int *restrict num_pop_particles
)
{
    for (int read = 0, write = 0, pop_idx = 0; read < *num_particles; read++)
    {
        if (pop_idx < *num_pop_particles && pop_indices[pop_idx] == read)
        {
            pop_idx++;
        }
        else
        {
            particle_indices[write] = particle_indices[read];
            write++;
        }
    }

    *num_particles -= *num_pop_particles;
    *num_pop_particles = 0;
}

/**
 * \brief Compute the acceleration between the particles 
 *        in the leaf and other nodes / leaves
 * 
 * \param a Array of acceleration vectors to be modified
 * \param x Array of position vectors
 * \param m Array of masses
 * \param G Gravitational constant
 * \param softening_length Softening length
 * \param opening_angle Opening angle
 * \param given_leaf Given leaf
 * \param root Root node of the octree
 * \param branch_record_given_leaf Branch record of the given leaf
 * \param branch_record_given_leaf_count Size of branch record
 * \param max_depth Maximum depth of the tree
 */
IN_FILE int _compute_acc_tree_walk(
    real *restrict a,
    const real *restrict x,
    const real *restrict m,
    const real G,
    const real softening_length,
    const real opening_angle,
    BarnesHutTreeLeaf *restrict given_leaf,
    BarnesHutTreeNode *restrict root,
    const int *restrict branch_record_given_leaf,
    const int branch_record_given_leaf_count,
    const int max_depth
)
{
    typedef struct BarnesHutAccTreeWalkStack
    {
        BarnesHutTreeNode *node;
        struct BarnesHutAccTreeWalkStack *last;
        int particle_indices[8 * MAX_NUM_PARTICLES_PER_LEAF];
        int num_particles[8];
        int processed_region;
    } BarnesHutAccTreeWalkStack;

    int return_code;

    BarnesHutTreeNode *node = root;

    /* Allocate a stack pool */
    const int stack_pool_size = max_depth;
    BarnesHutAccTreeWalkStack *stack_pool = malloc(stack_pool_size * sizeof(BarnesHutAccTreeWalkStack));
    if (!stack_pool)
    {
        return_code = ERROR_BARNES_HUT_COMPUTE_ACC_TREE_WALK_STACK_MEMORY_ALLOC;
        goto err_stack_pool_memory_alloc;
    }

    int stack_count = 1;
    BarnesHutAccTreeWalkStack *stack = &stack_pool[0];
    stack->node = root;
    stack->last = NULL;
    for (int i = 0; i < 8; i++)
    {
        stack->num_particles[i] = given_leaf->objects_count;
        for (int j = 0; j < given_leaf->objects_count; j++)
        {
            stack->particle_indices[i * MAX_NUM_PARTICLES_PER_LEAF + j] = given_leaf->indices[j];
        }
    }
    stack->processed_region = -1;

    int branch_record_count = 0;
    int *restrict branch_record = malloc(max_depth * sizeof(int));
    if (!branch_record)
    {
        return_code = ERROR_BARNES_HUT_COMPUTE_ACC_TREE_WALK_BRANCH_RECORD_MEMORY_ALLOC;
        goto err_branch_record_memory_alloc;
    }

    while (true)
    {
        for (int i = (stack->processed_region + 1); i < 8; i++)
        {
            BarnesHutTreeNode *child = node->children[i];
            BarnesHutTreeLeaf *leaf = node->leaves[i];

            // The node and leaf is empty
            if ((!child) && (!leaf))
            {
                stack->processed_region = i;
            }

            // Leaf exists
            else if (!child)
            {
                branch_record[branch_record_count] = i;
                branch_record_count++;

                bool is_same_leaf = _compare_branch_record(
                    branch_record_given_leaf,
                    branch_record_given_leaf_count,
                    branch_record,
                    branch_record_count
                );
                branch_record_count--;
                if (!is_same_leaf)
                {
                    _compute_acc_leaf_to_leaf(
                        a,
                        x,
                        m,
                        G,
                        softening_length,
                        &(stack->particle_indices[i * MAX_NUM_PARTICLES_PER_LEAF]),
                        stack->num_particles[i],
                        leaf
                    );
                }
                stack->processed_region = i;
            }

            // Node exists
            else
            {
                branch_record[branch_record_count] = i;
                branch_record_count++;

                bool is_included = _compare_branch_record(
                    branch_record_given_leaf,
                    branch_record_given_leaf_count,
                    branch_record,
                    branch_record_count
                );

                if (!is_included)
                {
                    real R[3];
                    const real COM_node[3] = {
                        child->center_of_mass[0],
                        child->center_of_mass[1],
                        child->center_of_mass[2]
                    };
                    const real m_node = child->total_mass;
                    const real box_width = child->box_width;

                    int pop_indices[MAX_NUM_PARTICLES_PER_LEAF];
                    int num_pop_particles = 0;
                    for (int j = 0; j < stack->num_particles[i]; j++)
                    {
                        const int idx_j = stack->particle_indices[i * MAX_NUM_PARTICLES_PER_LEAF + j];
                        
                        // Calculate \vec{R} and its norm
                        R[0] = x[idx_j * 3 + 0] - COM_node[0];
                        R[1] = x[idx_j * 3 + 1] - COM_node[1];
                        R[2] = x[idx_j * 3 + 2] - COM_node[2];

                        if (box_width / sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]) < opening_angle)
                        {
                            real temp_vec[3];
                            const real R_norm = sqrt(
                                R[0] * R[0] +
                                R[1] * R[1] +
                                R[2] * R[2] +
                                softening_length * softening_length
                            );
                            const real temp_value = G / (R_norm * R_norm * R_norm);
                            temp_vec[0] = temp_value * R[0];
                            temp_vec[1] = temp_value * R[1];
                            temp_vec[2] = temp_value * R[2];
                            a[idx_j * 3 + 0] -= temp_vec[0] * m_node;
                            a[idx_j * 3 + 1] -= temp_vec[1] * m_node;
                            a[idx_j * 3 + 2] -= temp_vec[2] * m_node;
                            pop_indices[num_pop_particles] = idx_j;
                            num_pop_particles++;
                        }
                    }

                    _pop_particles(
                        &(stack->particle_indices[i * MAX_NUM_PARTICLES_PER_LEAF]),
                        &(stack->num_particles[i]),
                        pop_indices,
                        &num_pop_particles
                    );

                    if (stack->num_particles[i] == 0)
                    {
                        stack->processed_region = i;
                        branch_record_count--;
                        break;
                    }
                }

                // Create a new stack item if the leaf is included
                if (stack_count >= stack_pool_size)
                {
                    return_code = ERROR_BARNES_HUT_COMPUTE_ACC_TREE_WALK_STACK_FULL;
                    goto err_stack_pool_full;
                }
                BarnesHutAccTreeWalkStack *new_item = &stack_pool[stack_count];
                stack_count++;

                new_item->node = child;
                new_item->last = stack;
                for (int j = 0; j < 8; j++)
                {
                    new_item->num_particles[j] = stack->num_particles[i];
                    for (int k = 0; k < stack->num_particles[i]; k++)
                    {
                        new_item->particle_indices[j * MAX_NUM_PARTICLES_PER_LEAF + k] = stack->particle_indices[i * MAX_NUM_PARTICLES_PER_LEAF + k];
                    }
                }
                new_item->processed_region = -1;

                stack->num_particles[i] = 0;    // Should be completely processed in the deeper level
                stack = new_item;
                node = child;

                break;
            }
        }

        if (stack->processed_region >= 7)
        {
            BarnesHutAccTreeWalkStack *parent_stack = stack->last;

            // Root node has no parent
            if (!parent_stack)
            {
                break;
            }

            stack = parent_stack;
            node = parent_stack->node;
            stack_count--;
            branch_record_count--;

            parent_stack->processed_region += 1;
        }
    }
    /* Free the memory */
    free(branch_record);
    free(stack_pool);

    return SUCCESS;

err_stack_pool_full:
err_branch_record_memory_alloc:
    free(branch_record);
err_stack_pool_memory_alloc:
    free(stack_pool);
    return return_code;
}

/**
 * \brief Compute the Barnes-Hut acceleration given the tree
 * 
 * \param a Array of acceleration vectors to be modified
 * \param x Array of position vectors
 * \param m Array of masses
 * \param G Gravitational constant
 * \param softening_length Softening length
 * \param opening_angle Opening angle
 * \param max_depth Maximum depth of the tree
 * \param root Root node of the octree
 * 
 * \retval SUCCESS If the computation is successful
 * \retval ERROR_BARNES_HUT_ACCELERATION_STEP_ACC_STACK_MEMORY_ALLOC 
 *         If memory allocation failed for the acceleration stack
 * \retval ERROR_BARNES_HUT_ACCELERATION_STEP_ACC_STACK_FULL 
 *         If the acceleration stack is full
 * \retval ERROR_BARNES_HUT_ACCELERATION_STEP_BRANCH_RECORD_MEMORY_ALLOC
 *         If memory allocation failed for the branch record
 * \retval error code if other errors occurred
 */
IN_FILE int _compute_acceleration(
    real *restrict a,
    const real *restrict x,
    const real *restrict m,
    const real G,
    const real softening_length,
    const real opening_angle,
    const int max_depth,
    BarnesHutTreeNode *root
)
{
    int return_code;

    BarnesHutTreeNode *node = root;
    
    /* Allocate a stack pool */
    const int stack_pool_size = max_depth;
    BarnesHutAccStack *stack_pool = malloc(stack_pool_size * sizeof(BarnesHutAccStack));
    if (!stack_pool)
    {
        return_code = ERROR_BARNES_HUT_COMPUTE_ACCELERATION_STACK_POOL_MEMORY_ALLOC;
        goto err_stack_pool_memory_alloc;
    }

    int stack_count = 1;
    BarnesHutAccStack *stack = &(stack_pool[0]);
    stack->node = root;
    stack->last = NULL;
    stack->processed_region = -1;

    // Keep track of the nodes that is in the same branch as the current node
    int branch_record_count = 0;
    int *restrict branch_record = malloc(max_depth * sizeof(int));
    if (!branch_record)
    {
        return_code = ERROR_BARNES_HUT_COMPUTE_ACCELERATION_BRANCH_RECORD_MEMORY_ALLOC;
        goto err_branch_record_memory_alloc;
    }

    while (true)
    {     
        /*
        *   First of all, we attempt to find a leaf node. Then, 
        *   (1) we calculate the acceleration between all 
        *       particles within the leaf. 
        *   (2) we calculate the acceleration between 
        *       the particles in the leaf and other nodes that 
        *       satisfy the condition s / d < opening_angle, 
        *       where s is the width of the node and d is the 
        *       distance between the particle and center of mass
        *       of the other node. In addition, we must note that the
        *       other node cannot be predecessor of the leaf, otherwise
        *       we are including the gravitational effect of the 
        *       particle due to itself, which is incorrect.
        */
        for (int i = (stack->processed_region + 1); i < 8; i++)
        {
            BarnesHutTreeNode *child = node->children[i];
            BarnesHutTreeLeaf *leaf = node->leaves[i];
            
            // The node and leaf is empty
            if ((!child) && (!leaf))
            {
                stack->processed_region = i;
            }

            // Leaf exists
            else if (!child)
            {
                branch_record[branch_record_count] = i;
                branch_record_count++;

                _compute_acc_single_leaf(
                    a,
                    x,
                    m,
                    G,
                    softening_length,
                    leaf
                );
                return_code = _compute_acc_tree_walk(
                    a,
                    x,
                    m,
                    G,
                    softening_length,
                    opening_angle,
                    leaf,
                    root,
                    branch_record,
                    branch_record_count,
                    max_depth
                );
                if (return_code != SUCCESS)
                {
                    goto err_compute_acc_tree_walk;
                }
                branch_record_count--;
                stack->processed_region = i;
            }

            // Node exists
            else
            {
                branch_record[branch_record_count] = i;
                branch_record_count++;

                // Create a new stack item
                if (stack_count >= stack_pool_size)
                {
                    return_code = ERROR_BARNES_HUT_COMPUTE_ACCELERATION_STACK_POOL_FULL;
                    goto err_stack_pool_full;
                }
                BarnesHutAccStack *new_item = &stack_pool[stack_count];
                stack_count++;

                new_item->node = child;
                new_item->last = stack;
                new_item->processed_region = -1;

                stack = new_item;
                node = child;

                break;
            }
        }

        if (stack->processed_region >= 7)
        {
            BarnesHutAccStack *parent_stack = stack->last;

            // Root node has no parent
            if (!parent_stack)
            {
                break;
            }

            parent_stack->processed_region += 1;
            stack = parent_stack;
            node = parent_stack->node;
            stack_count--;
            branch_record_count--;
        }
    }

    // Free the memory
    free(branch_record);
    free(stack_pool);
    return SUCCESS;

err_compute_acc_tree_walk:
err_stack_pool_full:
err_branch_record_memory_alloc:
    free(branch_record);
err_stack_pool_memory_alloc:
    free(stack_pool);
    return return_code;
}

WIN32DLL_API int acceleration_barnes_hut(
    real *restrict a,
    const System *restrict system,
    AccelerationParam *restrict acceleration_param
)
{
    int return_code;

    const int objects_count = system->objects_count;
    const real *restrict x = system->x;
    const real *restrict m = system->m;
    const real G = system->G;
    const real softening_length = acceleration_param->softening_length; 
    const real opening_angle = acceleration_param->opening_angle;

    /* Empty the input array */
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Find the width and center of the bounding box */
    real center[3];
    real width;
    _calculate_bounding_box(objects_count, x, center, &width);

    /* Construct the octree */
    // Allocate memory

    // root node
    BarnesHutTreeNode *restrict root = malloc(sizeof(BarnesHutTreeNode));
    if (root == NULL)
    {
        return_code = ERROR_BARNES_HUT_ROOT_MEMORY_ALLOC;
        goto err_root_memory;
    }
    _initialize_node(root, center, 0.0, width);

    // leaf pool
    BarnesHutTreeLeafPool *leaf_pool = malloc(sizeof(BarnesHutTreeNodePool));
    leaf_pool->pool_size = objects_count; // This value may not be optimal
    leaf_pool->leaves = NULL;
    leaf_pool->next = NULL;
    if (!leaf_pool)
    {
        return_code = ERROR_BARNES_HUT_LEAF_POOL_PTR_MEMORY_ALLOC;
        goto err_leaf_pool_ptr_memory;
    }
    leaf_pool->leaves = malloc(leaf_pool->pool_size * sizeof(BarnesHutTreeNode));
    if (!(leaf_pool->leaves))
    {
        return_code = ERROR_BARNES_HUT_LEAF_POOL_MEMORY_ALLOC;
        goto err_leaf_pool_memory;
    }

    // node pool
    BarnesHutTreeNodePool *node_pool = malloc(sizeof(BarnesHutTreeNodePool));
    node_pool->pool_size = objects_count; // This value may not be optimal
    node_pool->nodes = NULL;
    node_pool->next = NULL;
    if (!node_pool)
    {
        return_code = ERROR_BARNES_HUT_NODE_POOL_PTR_MEMORY_ALLOC;
        goto err_node_pool_ptr_memory;
    }
    node_pool->nodes = malloc(node_pool->pool_size * sizeof(BarnesHutTreeNode));
    if (!(node_pool->nodes))
    {
        return_code = ERROR_BARNES_HUT_NODE_POOL_MEMORY_ALLOC;
        goto err_node_pool_memory;
    }

    // Construct the octree
    int max_depth = 0;
    return_code = _construct_octree(
        objects_count,
        x,
        m,
        leaf_pool,
        node_pool,
        &max_depth,
        root
    );
    if (return_code != SUCCESS)
    {
        goto err_octree;
    }

    /* Calculate the center of mass */
    return_code = _compute_center_of_mass(x, m, max_depth, root);
    if (return_code != SUCCESS)
    {
        goto err_center_of_mass;
    }

    /* Calculate the acceleration */
    return_code = _compute_acceleration(a, x, m, G, softening_length, opening_angle, max_depth, root);
    if (return_code != SUCCESS)
    {
        goto err_acceleration;
    }

    /* Free the memory */
    while (node_pool != NULL)
    {
        BarnesHutTreeNodePool *next = node_pool->next;
        free(node_pool->nodes);
        free(node_pool);
        node_pool = next;
    }
    while (leaf_pool != NULL)
    {
        BarnesHutTreeLeafPool *next = leaf_pool->next;
        free(leaf_pool->leaves);
        free(leaf_pool);
        leaf_pool = next;
    }
    free(root);

    return SUCCESS;

err_acceleration:
err_center_of_mass:
err_octree:
err_node_pool_memory:
err_node_pool_ptr_memory:
    while (node_pool != NULL)
    {
        BarnesHutTreeNodePool *next = node_pool->next;
        free(node_pool->nodes);
        free(node_pool);
        node_pool = next;
    }
err_leaf_pool_memory:
err_leaf_pool_ptr_memory:
    while (leaf_pool != NULL)
    {
        BarnesHutTreeLeafPool *next = leaf_pool->next;
        free(leaf_pool->leaves);
        free(leaf_pool);
        leaf_pool = next;
    }
err_root_memory:
    free(root);
    
    return return_code;
}
