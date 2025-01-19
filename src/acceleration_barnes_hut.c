#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "acceleration_barnes_hut.h"

// For debug
// void visualize_octree_nodes(BarnesHutTreeNode *node, int depth) 
// {
//     if (node == NULL) 
//     {
//         return;
//     }

//     // Indentation for better visualization
//     for (int i = 0; i < depth; i++) 
//     {
//         printf("  ");
//     }

//     // Print node information
//     printf("Node (depth %d): is_leaf=%s, index=%d, mass=%.10f, center=(%.10f, %.10f, %.10f), width=%.2f, intern_obj_count=%d\n",
//            depth, node->index >= 0 ? "true" : "false", node->index, node->total_mass,
//            node->center_of_mass[0], node->center_of_mass[1], node->center_of_mass[2],
//            node->box_width, node->internal_nodes_objects_count);

//     // Recursively visualize children
//     for (int i = 0; i < 8; i++) 
//     {
//         visualize_octree_nodes(node->children[i], depth + 1);
//     }
// }

WIN32DLL_API void acceleration_barnes_hut(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length,
    real barnes_hut_theta
)
{
    // Empty the input array
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
    // root node
    BarnesHutTreeNode *restrict root = malloc(sizeof(BarnesHutTreeNode));
    if (root == NULL) 
    {
        fprintf(stderr, "Error: Failed to allocate memory for the root node in acceleration_barnes_hut.\n");
        goto err_root_memory;
    }
    root->index = -1;
    root->total_mass = 0.0;
    root->center_of_mass[0] = center[0];
    root->center_of_mass[1] = center[1];
    root->center_of_mass[2] = center[2];
    root->box_width = width;
    root->internal_nodes_objects_count = 0;
    root->children[0] = NULL;
    root->children[1] = NULL;
    root->children[2] = NULL;
    root->children[3] = NULL;
    root->children[4] = NULL;
    root->children[5] = NULL;
    root->children[6] = NULL;
    root->children[7] = NULL;

    // leaf node pool
    BarnesHutTreeNodePool *leaf_node_pool = malloc(sizeof(BarnesHutTreeNodePool));
    if (leaf_node_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for leaf nodes in acceleration_barnes_hut.\n");
        goto err_leaf_pool_memory;
    }
    leaf_node_pool->pool_size = objects_count;
    leaf_node_pool->node_pool = malloc(leaf_node_pool->pool_size * sizeof(BarnesHutTreeNode));
    leaf_node_pool->next = NULL;
    if (leaf_node_pool->node_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for leaf nodes in acceleration_barnes_hut.\n");
        goto err_leaf_pool_memory;
    }

    // internal node pool
    BarnesHutTreeNodePool *internal_node_pool = malloc(sizeof(BarnesHutTreeNodePool));
    if (internal_node_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for internal nodes in acceleration_barnes_hut.\n");
        goto err_internal_pool_memory;
    }
    internal_node_pool->pool_size = leaf_node_pool->pool_size;
    internal_node_pool->node_pool = malloc(internal_node_pool->pool_size * sizeof(BarnesHutTreeNode));
    internal_node_pool->next = NULL;
    if (internal_node_pool->node_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for internal nodes in acceleration_barnes_hut.\n");
        goto err_internal_pool_memory;
    }
    
    // Construct the octree
    int actual_interval_nodes_count = 0;
    if (_barnes_hut_construct_octree(objects_count, x, m, width, leaf_node_pool, internal_node_pool, &actual_interval_nodes_count, root) == 1)
    {
        goto err_octree_memory;
    }

    // Shorten the tree if internal node have less than 8 leaves
    if (_barnes_hut_shorten_tree(&actual_interval_nodes_count, root) == 1)
    {
        goto err_octree_memory;
    }

    /* Calculate the center of mass */
    if (_barnes_hut_compute_center_of_mass(actual_interval_nodes_count, root) == 1)
    {
        goto err_center_of_mass_memory;
    }

    /* Calculate the acceleration */
    if (_barnes_hut_acceleration(a, G, softening_length, barnes_hut_theta, actual_interval_nodes_count, root) == 1)
    {
        goto err_acceleration_memory;
    }

    // Free the memory
    free(leaf_node_pool->node_pool);
    free(leaf_node_pool);
    while (internal_node_pool != NULL)
    {
        BarnesHutTreeNodePool *next = internal_node_pool->next;
        free(internal_node_pool->node_pool);
        free(internal_node_pool);
        internal_node_pool = next;
    }
    free(root);
    
    return;

err_acceleration_memory:
err_center_of_mass_memory:
err_octree_memory:
err_internal_pool_memory:
    free(internal_node_pool);
err_leaf_pool_memory:
    free(leaf_node_pool);
err_root_memory:
    free(root);
    return;
//    return 1;
}

WIN32DLL_API void _calculate_bounding_box(
    int objects_count,
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

WIN32DLL_API int _barnes_hut_check_region(
    real x,
    real y,
    real z,
    real center_x,
    real center_y,
    real center_z
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

WIN32DLL_API int _barnes_hut_construct_octree(
    int objects_count,
    const real *restrict x,
    const real *restrict m,
    real width,
    BarnesHutTreeNodePool *leaf_node_pool,
    BarnesHutTreeNodePool *internal_node_pool,
    int *restrict actual_interval_nodes_count,
    BarnesHutTreeNode *restrict root
)
{
    // Note: This variable is used to keep track of whether
    //       the node pool is full. It will be reset to 0 
    //       when the node pool is full and expanded. Thus,
    //       don't trust the value of this variable.
    int current_internal_nodes_count = 0;
    for (int i = 0; i < objects_count; i++)
    {
        int depth = 1;
        BarnesHutTreeNode *current_node = root;
        BarnesHutTreeNode *leaf_node = &(leaf_node_pool->node_pool[i]);

        leaf_node->index = i;
        leaf_node->total_mass = m[i];
        leaf_node->center_of_mass[0] = x[i * 3 + 0];
        leaf_node->center_of_mass[1] = x[i * 3 + 1];
        leaf_node->center_of_mass[2] = x[i * 3 + 2];
        // Technically not used
        // leaf_node->box_width = 0.0;
        // leaf_node->internal_nodes_objects_count = 0;
        // leaf_node->children[0] = NULL;
        // leaf_node->children[1] = NULL;
        // leaf_node->children[2] = NULL;
        // leaf_node->children[3] = NULL;
        // leaf_node->children[4] = NULL;
        // leaf_node->children[5] = NULL;
        // leaf_node->children[6] = NULL;
        // leaf_node->children[7] = NULL;

        BarnesHutTreeNode *collider_leaf = NULL;
        
        while (true)
        {
            /*
            *   Note: In the following parameters,
            *         the center_of_mass is just
            *         the center of the node, not
            *         the center of mass. It will
            *         be replaced in the next step.
            */
            int region = _barnes_hut_check_region(
                x[i * 3 + 0],
                x[i * 3 + 1],
                x[i * 3 + 2],
                current_node->center_of_mass[0],
                current_node->center_of_mass[1],
                current_node->center_of_mass[2]
            );

            (current_node->internal_nodes_objects_count)++;
            
            if (collider_leaf == NULL)
            {
                // The node is empty
                if (current_node->children[region] == NULL)
                {
                    current_node->children[region] = leaf_node;
                    break;
                }

                // The node is already occupied by a leaf. (>= 0 means it is a leaf)
                else if ((current_node->children[region]->index) >= 0)
                {
                    collider_leaf = current_node->children[region];

                    // Create a new internal node
                    if (current_internal_nodes_count >= internal_node_pool->pool_size)
                    {
                        internal_node_pool->next = malloc(sizeof(BarnesHutTreeNodePool));
                        if (internal_node_pool->next == NULL)
                        {
                            fprintf(stderr, "Error: Failed to allocate memory for internal nodes in _barnes_hut_construct_octree.\n");
                            goto err_memory;
                        }
                        internal_node_pool->next->pool_size = internal_node_pool -> pool_size * 2;
                        internal_node_pool->next->node_pool = malloc(internal_node_pool->next->pool_size * sizeof(BarnesHutTreeNode));
                        if (internal_node_pool->next->node_pool == NULL)
                        {
                            fprintf(stderr, "Error: Failed to allocate memory for internal nodes in _barnes_hut_construct_octree.\n");
                            goto err_memory;
                        }
                        internal_node_pool->next->next = NULL;
                        internal_node_pool = internal_node_pool->next;

                        current_internal_nodes_count = 0;
                    }
                    BarnesHutTreeNode *new_node = &(internal_node_pool->node_pool[current_internal_nodes_count]);
                    (*actual_interval_nodes_count)++;
                    current_internal_nodes_count++;

                    new_node->index = -1;
                    new_node->total_mass = collider_leaf->total_mass + m[i];
                    new_node->box_width = width / ((real) fast_pow_of_2(depth));
                    new_node->internal_nodes_objects_count = 1;
                    new_node->children[0] = NULL;
                    new_node->children[1] = NULL;
                    new_node->children[2] = NULL;
                    new_node->children[3] = NULL;
                    new_node->children[4] = NULL;
                    new_node->children[5] = NULL;
                    new_node->children[6] = NULL;
                    new_node->children[7] = NULL;

                    // Calculate the center of the node (not center of mass) when constructing the octree
                    // It will be replaced with the center of mass in the next step
                    if (region & 1)
                    {
                        new_node->center_of_mass[0] = current_node->center_of_mass[0] + width / ((real) fast_pow_of_2(depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[0] = current_node->center_of_mass[0] - width / ((real) fast_pow_of_2(depth + 1));
                    }

                    if (region & 2)
                    {
                        new_node->center_of_mass[1] = current_node->center_of_mass[1] + width / ((real) fast_pow_of_2(depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[1] = current_node->center_of_mass[1] - width / ((real) fast_pow_of_2(depth + 1));
                    }

                    if (region & 4)
                    {
                        new_node->center_of_mass[2] = current_node->center_of_mass[2] + width / ((real) fast_pow_of_2(depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[2] = current_node->center_of_mass[2] - width / ((real) fast_pow_of_2(depth + 1));
                    }

                    current_node->children[region] = new_node;
                    current_node = new_node;
                    depth++;
                }

                // The node is already occupied by a internal node
                else
                {   
                    current_node = current_node->children[region];
                    current_node->total_mass += m[i];
                    depth++;
                }
            }
            else
            {
                int collider_region = _barnes_hut_check_region(
                    collider_leaf->center_of_mass[0],
                    collider_leaf->center_of_mass[1],
                    collider_leaf->center_of_mass[2],
                    current_node->center_of_mass[0],
                    current_node->center_of_mass[1],
                    current_node->center_of_mass[2]
                );

                if (region != collider_region)
                {
                    current_node->children[region] = leaf_node;
                    current_node->children[collider_region] = collider_leaf;
                    break;
                }

                // Collide again
                else
                {
                    // Create a new internal node
                    if (current_internal_nodes_count >= internal_node_pool->pool_size)
                    {
                        internal_node_pool->next = malloc(sizeof(BarnesHutTreeNodePool));
                        if (internal_node_pool->next == NULL)
                        {
                            fprintf(stderr, "Error: Failed to allocate memory for internal nodes in _barnes_hut_construct_octree.\n");
                            goto err_memory;
                        }
                        internal_node_pool->next->pool_size = internal_node_pool -> pool_size * 2;
                        internal_node_pool->next->node_pool = malloc(internal_node_pool->next->pool_size * sizeof(BarnesHutTreeNode));
                        if (internal_node_pool->next->node_pool == NULL)
                        {
                            fprintf(stderr, "Error: Failed to allocate memory for internal nodes in _barnes_hut_construct_octree.\n");
                            goto err_memory;
                        }
                        internal_node_pool->next->next = NULL;
                        internal_node_pool = internal_node_pool->next;

                        current_internal_nodes_count = 0;
                    }
                    BarnesHutTreeNode *new_node = &(internal_node_pool->node_pool[current_internal_nodes_count]);
                    (*actual_interval_nodes_count)++;
                    current_internal_nodes_count++;

                    new_node->index = -1;
                    new_node->total_mass = current_node->total_mass;
                    new_node->center_of_mass[0] = current_node->center_of_mass[0];
                    new_node->center_of_mass[1] = current_node->center_of_mass[1];
                    new_node->center_of_mass[2] = current_node->center_of_mass[2];
                    new_node->box_width = width / ((real) fast_pow_of_2(depth));
                    new_node->internal_nodes_objects_count = 1;
                    new_node->children[0] = NULL;
                    new_node->children[1] = NULL;
                    new_node->children[2] = NULL;
                    new_node->children[3] = NULL;
                    new_node->children[4] = NULL;
                    new_node->children[5] = NULL;
                    new_node->children[6] = NULL;
                    new_node->children[7] = NULL;

                    // Calculate the center of the node when constructing the octree
                    // It will be replaced with the center of mass later
                    if (region & 1)
                    {
                        new_node->center_of_mass[0] = current_node->center_of_mass[0] + width / ((real) fast_pow_of_2(depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[0] = current_node->center_of_mass[0] - width / ((real) fast_pow_of_2(depth + 1));
                    }

                    if (region & 2)
                    {
                        new_node->center_of_mass[1] = current_node->center_of_mass[1] + width / ((real) fast_pow_of_2(depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[1] = current_node->center_of_mass[1] - width / ((real) fast_pow_of_2(depth + 1));
                    }
                    
                    if (region & 4)
                    {
                        new_node->center_of_mass[2] = current_node->center_of_mass[2] + width / ((real) fast_pow_of_2(depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[2] = current_node->center_of_mass[2] - width / ((real) fast_pow_of_2(depth + 1));
                    }

                    current_node->children[region] = new_node;
                    current_node = new_node;
                    depth++;
                }
            }
        }
    }

    return 0;

err_memory:
    return 1;
}

WIN32DLL_API int _barnes_hut_shorten_tree(
    int *restrict actual_interval_nodes_count,
    BarnesHutTreeNode *restrict root
)
{
    if (root->internal_nodes_objects_count <= 8)
    {
        return 0;
    }

    typedef struct BarnesHutShortenTreeStack
    {
        BarnesHutTreeNode *node;
        struct BarnesHutShortenTreeStack *last;
        int processed_region;
    } BarnesHutShortenTreeStack;

    BarnesHutTreeNode *current_node = root;

    int shorten_stack_pool_size = *actual_interval_nodes_count + root->internal_nodes_objects_count + 1;
    BarnesHutShortenTreeStack *shorten_tree_stack_pool = malloc(shorten_stack_pool_size * sizeof(BarnesHutShortenTreeStack));
    if (shorten_tree_stack_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for stack in _barnes_hut_shorten_tree.\n");
        goto err_memory;
    }
    int current_stack_count = 1;
    BarnesHutShortenTreeStack *stack = &shorten_tree_stack_pool[0];

    stack->node = root;
    stack->last = NULL;
    stack->processed_region = -1;

    bool is_shorten = false;
    BarnesHutTreeNode *shorten_leaves[8];
    int num_leaves_to_find;
    int found_leaves_count;
    while (true)
    {
        for (int i = (stack->processed_region + 1); i < 8; i++)
        {
            BarnesHutTreeNode *child_i = current_node->children[i];

            // Check if the node is empty 
            // or occupied by a leaf (index >= 0 means it is a leaf) 
            if (
                (child_i == NULL)
                || (child_i->index >= 0)
            )
            {
                stack->processed_region = i;
                continue;
            }

            else
            {
                // Create a new stack item
                if (current_stack_count >= shorten_stack_pool_size)
                {
                    fprintf(stderr, "Error: The stack pool is full in _barnes_hut_shorten_tree.\n");
                    goto err_memory;
                }
                BarnesHutShortenTreeStack *new_item = &shorten_tree_stack_pool[current_stack_count];
                current_stack_count++;

                new_item->node = child_i;
                new_item->last = stack;
                new_item->processed_region = -1;

                stack = new_item;
                current_node = child_i;
                break;
            }
        }

        if (stack->processed_region >= 7)
        {
            // Root node has no parent
            if (stack->last == NULL)
            {
                break;
            }

            stack = stack->last;
            (stack->processed_region)++;
            current_node = stack->node;
        }

        else if (current_node->internal_nodes_objects_count <= 8)
        {
            for (int i = 0; i < 8; i++)
            {
                BarnesHutTreeNode *temp_child = current_node->children[i];
                // Check if there is internal node
                if ((temp_child != NULL) && (temp_child->index < 0))
                {
                    is_shorten = true;
                    break;
                }
            }
        }

        if (!is_shorten)
        {
            continue;
        }

        // Find the leaves
        num_leaves_to_find = current_node->internal_nodes_objects_count;
        found_leaves_count = 0;
        BarnesHutShortenTreeStack *find_leaves_stack = stack;
        BarnesHutTreeNode *find_leaves_current_node = current_node;
        while (true)
        {
            for (int i = find_leaves_stack->processed_region + 1; i < 8; i++)
            {
                BarnesHutTreeNode *temp_child = find_leaves_current_node->children[i];
                
                // Check if the node is empty
                if (temp_child == NULL)
                {
                    find_leaves_stack->processed_region = i;
                    continue;
                }

                // Check if the node is occupied by a leaf
                else if (temp_child->index >= 0)
                {
                    shorten_leaves[found_leaves_count] = temp_child;
                    found_leaves_count++;
                    find_leaves_stack->processed_region = i;
                    break;
                }

                else
                {
                    // Create a new stack item
                    if (current_stack_count >= shorten_stack_pool_size)
                    {
                        fprintf(stderr, "Error: The stack pool is full in _barnes_hut_shorten_tree2.\n");
                        goto err_memory;
                    }
                    BarnesHutShortenTreeStack *new_item = &shorten_tree_stack_pool[current_stack_count];
                    current_stack_count++;

                    new_item->node = temp_child;
                    new_item->last = find_leaves_stack;
                    new_item->processed_region = -1;

                    find_leaves_stack = new_item;
                    find_leaves_current_node = temp_child;
                    break;
                }
            }

            if (found_leaves_count >= num_leaves_to_find)
            {
                break;
            }

            if (find_leaves_stack->processed_region >= 7)
            {
                BarnesHutShortenTreeStack *parent = find_leaves_stack->last;
                if (parent == NULL)
                {
                    break;
                }

                find_leaves_stack = parent;
                find_leaves_current_node = parent->node;
                (find_leaves_stack->processed_region)++;
            }
        }

        if (found_leaves_count < num_leaves_to_find)
        {
            fprintf(stderr, "Error: Failed to find all leaves in _barnes_hut_shorten_tree.\n");
            goto err_shorten_tree;
        }

        for (int i = 0, placed_leaves_count = 0; (i < 8) && (placed_leaves_count < found_leaves_count); i++)
        {
            current_node->children[i] = shorten_leaves[placed_leaves_count];
            placed_leaves_count++;
        }
        for (int i = num_leaves_to_find; i < 8; i++)
        {
            current_node->children[i] = NULL;
        }

        is_shorten = false;
        stack->processed_region = 7;
    }

    // Free the stack pool
    free(shorten_tree_stack_pool);

    return 0;

err_shorten_tree:
    free(shorten_tree_stack_pool);
err_memory:
    return 1;
}

WIN32DLL_API int _barnes_hut_compute_center_of_mass(
    int actual_interval_nodes_count,
    BarnesHutTreeNode *restrict root
)
{
    typedef struct BarnesHutCOMStack
    {
        BarnesHutTreeNode *node;
        struct BarnesHutCOMStack *last;
        real sum_of_mass_times_distance[3];
        int processed_region;
    } BarnesHutCOMStack;

    BarnesHutTreeNode *current_node = root;

    int COM_stack_pool_size = actual_interval_nodes_count + 1;
    BarnesHutCOMStack *COM_stack_pool = malloc(COM_stack_pool_size * sizeof(BarnesHutCOMStack));
    if (COM_stack_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for stack in _barnes_hut_compute_center_of_mass.\n");
        goto err_memory;
    }
    int current_stack_count = 1;
    BarnesHutCOMStack *stack = &COM_stack_pool[0];

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
            BarnesHutTreeNode *child_i = current_node->children[i];

            // The node is empty
            if (child_i == NULL)
            {
                stack->processed_region = i;
                continue;
            }

            // The node is occupied by a leaf (>= 0 means it is a leaf)
            else if (child_i->index >= 0)
            {
                real child_mass = child_i->total_mass;
                stack->sum_of_mass_times_distance[0] += child_mass * child_i->center_of_mass[0];
                stack->sum_of_mass_times_distance[1] += child_mass * child_i->center_of_mass[1];
                stack->sum_of_mass_times_distance[2] += child_mass * child_i->center_of_mass[2];

                stack->processed_region = i;
                continue;
            }
            else
            {
                // Create a new stack item
                if (current_stack_count >= COM_stack_pool_size)
                {
                    fprintf(stderr, "Error: The stack pool is full in _barnes_hut_compute_center_of_mass.\n");
                    goto err_memory;
                }
                BarnesHutCOMStack *new_item = &COM_stack_pool[current_stack_count];
                current_stack_count++;

                new_item->node = child_i;
                new_item->last = stack;
                new_item->sum_of_mass_times_distance[0] = 0.0;
                new_item->sum_of_mass_times_distance[1] = 0.0;
                new_item->sum_of_mass_times_distance[2] = 0.0;
                new_item->processed_region = -1;

                stack = new_item;
                current_node = child_i;
                break;
            }
        }
    
        if (stack->processed_region >= 7)
        {
            BarnesHutCOMStack *parent = stack->last;
            
            // Root node has no parent
            if (parent == NULL)
            {
                break;
            }

            real total_mass = current_node->total_mass;

            current_node->center_of_mass[0] = stack->sum_of_mass_times_distance[0] / total_mass;
            current_node->center_of_mass[1] = stack->sum_of_mass_times_distance[1] / total_mass;
            current_node->center_of_mass[2] = stack->sum_of_mass_times_distance[2] / total_mass;

            BarnesHutTreeNode *parent_node = parent->node;
            parent->sum_of_mass_times_distance[0] += total_mass * current_node->center_of_mass[0];
            parent->sum_of_mass_times_distance[1] += total_mass * current_node->center_of_mass[1];
            parent->sum_of_mass_times_distance[2] += total_mass * current_node->center_of_mass[2];
            parent->processed_region++;

            stack = parent;
            current_node = parent_node;
        }
    }

    // Free the stack pool
    free(COM_stack_pool);

    return 0;

err_memory:
    return 1;
}

WIN32DLL_API int _barnes_hut_acceleration(
    real *restrict a,
    real G,
    real softening_length,
    real theta,
    int actual_interval_nodes_count,
    BarnesHutTreeNode *restrict root
)
{
    BarnesHutTreeNode *current_acc_node = root;

    typedef struct BarnesHutAccStack
    {
        BarnesHutTreeNode *node;
        struct BarnesHutAccStack *last;
        int processed_region;
    } BarnesHutAccStack;

    typedef struct BarnesHutSameBranchNode
    {
        BarnesHutTreeNode *node;
        struct BarnesHutSameBranchNode *next;
    } BarnesHutSameBranchNode;
    
    /* Allocate a stack pool for the outer loop */
    int acc_stack_pool_size = actual_interval_nodes_count + 1;
    BarnesHutAccStack *acc_stack_pool = malloc(acc_stack_pool_size * sizeof(BarnesHutAccStack));
    if (acc_stack_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for stack in _barnes_hut_acceleration.\n");
        goto err_acc_stack_pool_init_memory;
    }

    /* Allocate a stack pool for the inner loop */
    int obj_stack_pool_size = actual_interval_nodes_count + 1;
    BarnesHutAccStack *obj_stack_pool = malloc(obj_stack_pool_size * sizeof(BarnesHutAccStack));
    if (obj_stack_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for stack in _barnes_hut_acceleration.\n");
        goto err_obj_stack_pool_init_memory;
    }

    int current_acc_stack_count = 1;
    BarnesHutAccStack *acc_stack = &acc_stack_pool[0];

    acc_stack->node = current_acc_node;
    acc_stack->last = NULL;
    acc_stack->processed_region = -1;

    bool found_leaf;
    int acc_object_index = 0;
    BarnesHutTreeNode *current_acc_leaf;

    // Keep track of the nodes that is in the same branch as the current node
    int same_branch_node_pool_size = actual_interval_nodes_count + 1;
    BarnesHutSameBranchNode *acc_same_branch_node_pool = malloc(same_branch_node_pool_size * sizeof(BarnesHutSameBranchNode));
    if (acc_same_branch_node_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for acc_same_branch_nodes in _barnes_hut_acceleration.\n");
        goto err_same_branch_node_pool_init_memory;
    }
    int current_same_branch_nodes_count = 1;
    BarnesHutSameBranchNode *current_same_branch_node = &acc_same_branch_node_pool[0];
    current_same_branch_node->node = root;
    BarnesHutSameBranchNode *same_branch_node_root = current_same_branch_node;

    while (true)
    {     
        /*
        *   First of all, we attempt to find a leaf node.
        *   Then, we calculate the acceleration between 
        *   the leaf node and all other nodes that satisfy
        *   s / d < theta, where s is the width of the node
        *   and d is the distance between the center of mass
        *   of the node and the leaf node.
        */
        found_leaf = false;
        current_same_branch_nodes_count = 0;
        for (int i = (acc_stack->processed_region + 1); i < 8; i++)
        {   
            BarnesHutTreeNode *child_i = current_acc_node->children[i];

            // The node is empty
            if (child_i == NULL)
            {
                acc_stack->processed_region = i;
                continue;
            }

            // The node is occupied by a leaf (>= 0 means it is a leaf)
            else if (child_i->index >= 0)
            {
                current_acc_leaf = child_i;
                acc_object_index = current_acc_leaf->index;
                found_leaf = true;
                acc_stack->processed_region = i;

                if (current_same_branch_nodes_count >= same_branch_node_pool_size)
                {
                    fprintf(stderr, "Error: The same branch node pool is full in _barnes_hut_acceleration.\n");
                    goto err_memory;
                }
                current_same_branch_node->next = &acc_same_branch_node_pool[current_same_branch_nodes_count];
                current_same_branch_node = current_same_branch_node->next;
                current_same_branch_node->node = current_acc_leaf;
                current_same_branch_nodes_count++;
                break;
            }

            // The node is occupied by a internal node
            else
            {
                // Create a new stack item
                if (current_acc_stack_count >= acc_stack_pool_size)
                {
                    fprintf(stderr, "Error: The stack pool is full in _barnes_hut_acceleration.\n");
                    goto err_memory;
                }
                BarnesHutAccStack *new_item = &acc_stack_pool[current_acc_stack_count];
                current_acc_stack_count++;

                new_item->node = child_i;
                new_item->last = acc_stack;
                new_item->processed_region = -1;

                acc_stack = new_item;
                current_acc_node = child_i;

                if (current_same_branch_nodes_count >= same_branch_node_pool_size)
                {
                    fprintf(stderr, "Error: The same branch node pool is full in _barnes_hut_acceleration.\n");
                    goto err_memory;
                }
                current_same_branch_node->next = &acc_same_branch_node_pool[current_same_branch_nodes_count];
                current_same_branch_node = current_same_branch_node->next;
                current_same_branch_node->node = current_acc_node;
                current_same_branch_nodes_count++;
                break;
            }
        }
        
        if (found_leaf)
        {
            BarnesHutTreeNode *current_obj_node = root;

            // Note: This variable is used to keep track of whether
            //       the stack pool is full. It will be reset to 0
            //       when the stack is full and expanded. Thus, don't
            //       trust the value of this variable.
            int current_obj_stack_count = 0;

            BarnesHutAccStack *obj_stack = &obj_stack_pool[current_obj_stack_count];
            current_obj_stack_count++;

            obj_stack->node = root;
            obj_stack->last = NULL;
            obj_stack->processed_region = -1;

            // Calculate acceleration
            while (true)
            {
                for (int j = (obj_stack->processed_region + 1); j < 8; j++)
                {
                    BarnesHutTreeNode *child_j = current_obj_node->children[j];

                    // The node is empty
                    if (child_j == NULL)
                    {
                        obj_stack->processed_region = j;
                        continue;
                    }

                    // Check if the node is in the same branch as the current node
                    bool is_same_branch = false;
                    current_same_branch_node = same_branch_node_root;
                    for (int k = 0; k < current_same_branch_nodes_count; k++)
                    { 
                        if (current_same_branch_node->node == child_j)
                        {
                            is_same_branch = true;
                            break;
                        }
                        current_same_branch_node = current_same_branch_node->next;
                    }
                    if (is_same_branch)
                    {
                        obj_stack->processed_region = j;
                        continue;
                    }
                    
                    else
                    { 
                        real R[3];
                        real R_norm;

                        // Calculate \vec{R} and its norm
                        R[0] = current_acc_leaf->center_of_mass[0] - child_j->center_of_mass[0];
                        R[1] = current_acc_leaf->center_of_mass[1] - child_j->center_of_mass[1];
                        R[2] = current_acc_leaf->center_of_mass[2] - child_j->center_of_mass[2];

                        R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2] + softening_length * softening_length);

                        // Check if the node satisfies the condition (s / d < theta) or it is a leaf
                        if (((child_j->box_width / R_norm) < theta) || (child_j->index >= 0))
                        {
                            if (child_j->index != acc_object_index)
                            {
                                _barnes_hut_helper_acceleration_pair(current_acc_leaf, child_j, a, G, R, R_norm);
                            }
                            obj_stack->processed_region = j;
                            break;
                        }
                        else
                        {
                            // Create a new stack item
                            if (current_obj_stack_count >= (obj_stack_pool_size))
                            {
                                fprintf(stderr, "Error: The stack pool is full in _barnes_hut_acceleration.\n");
                                goto err_memory;
                            }

                            BarnesHutAccStack *new_item = &obj_stack_pool[current_obj_stack_count];
                            current_obj_stack_count++;

                            new_item->node = child_j;
                            new_item->last = obj_stack;
                            new_item->processed_region = -1;

                            obj_stack = new_item;
                            current_obj_node = child_j;
                            break;
                        }
                    }
                }
                
                if (obj_stack->processed_region >= 7)
                {
                    BarnesHutAccStack *parent = obj_stack->last;

                    // If the parent is NULL, it means that the current node is the root node
                    if (parent == NULL)
                    {
                        break;
                    }

                    parent->processed_region++;
                    obj_stack = parent;
                    current_obj_node = parent->node;
                }   
            }
        }

        if (acc_stack->processed_region >= 7)
        {
            BarnesHutAccStack *parent = acc_stack->last;

            // If the parent is NULL, it means that the current node is the root node
            if (parent == NULL)
            {
                break;
            }

            parent->processed_region++;
            acc_stack = parent;
            current_acc_node = parent->node;
        }
    }

    // Free the stack pool
    free(acc_stack_pool);
    free(obj_stack_pool);
    free(acc_same_branch_node_pool);    

    return 0;

err_memory:
err_same_branch_node_pool_init_memory:
    free(acc_same_branch_node_pool);
err_obj_stack_pool_init_memory:
    free(obj_stack_pool);
err_acc_stack_pool_init_memory:
    free(acc_stack_pool);
    return 1;
}

WIN32DLL_API void _barnes_hut_helper_acceleration_pair(
    BarnesHutTreeNode *restrict current_acc_leaf,
    BarnesHutTreeNode *restrict current_obj_leaf,
    real *restrict a,
    real G,
    real *restrict R,
    real R_norm
)
{
    real temp_value;

    // Calculate the acceleration
    int acc_object_index = current_acc_leaf->index;
    real object_2_mass = current_obj_leaf->total_mass;
    temp_value = G * object_2_mass / (R_norm * R_norm * R_norm);
    a[acc_object_index * 3 + 0] -= temp_value * R[0];
    a[acc_object_index * 3 + 1] -= temp_value * R[1];
    a[acc_object_index * 3 + 2] -= temp_value * R[2];
}