#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <time.h>

#include "common.h"
#include "acceleration.h"
#ifdef USE_CUDA
    #include "acceleration_cuda.cuh"
#endif

WIN32DLL_API int get_acceleration_method_flag(
    const char *restrict acceleration_method
)
{
    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        return ACCELERATION_METHOD_PAIRWISE;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        return ACCELERATION_METHOD_MASSLESS;
    }
    else if (strcmp(acceleration_method, "barnes-hut") == 0)
    {
        return ACCELERATION_METHOD_BARNES_HUT;
    }
    #ifdef USE_CUDA
        else if (strcmp(acceleration_method, "pairwise_cuda") == 0)
        {
            return ACCELERATION_METHOD_PAIRWISE_CUDA;
        }
        else if (strcmp(acceleration_method, "pairwise_float_cuda") == 0)
        {
            return ACCELERATION_METHOD_PAIRWISE_FLOAT_CUDA;
        }
        else if (strcmp(acceleration_method, "pairwise_float_comp_sum_cuda") == 0)
        {
            return ACCELERATION_METHOD_PAIRWISE_FLOAT_COMP_SUM_CUDA;
        }
    #endif
    else
    {
        return -1;
    }
}

WIN32DLL_API void acceleration(
    int acceleration_method_flag,
    int objects_count,
    real *restrict x,
    real *restrict v,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length,
    real barnes_hut_theta
)
{
    (void) v;   // To suppress the warning of unused variable
    switch (acceleration_method_flag)
    {
        // Pairwise acceleration
        case ACCELERATION_METHOD_PAIRWISE:
            acceleration_pairwise(objects_count, x, a, m, G, softening_length);
            break;
        // Massless acceleration
        case ACCELERATION_METHOD_MASSLESS:
            acceleration_massless(objects_count, x, a, m, G, softening_length);
            break;
        // Barnes-Hut acceleration
        case ACCELERATION_METHOD_BARNES_HUT:
            acceleration_barnes_hut(objects_count, x, a, m, G, softening_length, barnes_hut_theta);
            break;
        #ifdef USE_CUDA
            // CUDA Pairwise acceleration
            case ACCELERATION_METHOD_PAIRWISE_CUDA:
                acceleration_pairwise_cuda(objects_count, x, a, m, G, softening_length);
                break;
            // CUDA Pairwise acceleration with single precision
            case ACCELERATION_METHOD_PAIRWISE_FLOAT_CUDA:
                acceleration_pairwise_float_cuda(objects_count, x, a, m, G, softening_length);
                break;
        #endif
        default:
            fprintf(stderr, "Warning: Invalid acceleration method detected. Ignoring acceleration calculation.\n");
            break;
    }
}

WIN32DLL_API void acceleration_pairwise(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length
)
{   
    real R_norm;
    real temp_value;
    real temp_vec[3];
    real R[3];
    
    // Empty the input array
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    for(int i = 0; i < objects_count; i++)
    {
        for(int j = i + 1; j < objects_count; j++)
        {
            // Calculate \vec{R} and its norm
            R[0] = x[i * 3 + 0] - x[j * 3 + 0];
            R[1] = x[i * 3 + 1] - x[j * 3 + 1];
            R[2] = x[i * 3 + 2] - x[j * 3 + 2];
            R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2] + softening_length * softening_length);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[i * 3 + 0] -= temp_vec[0] * m[j];
            a[i * 3 + 1] -= temp_vec[1] * m[j];
            a[i * 3 + 2] -= temp_vec[2] * m[j];
            a[j * 3 + 0] += temp_vec[0] * m[i];
            a[j * 3 + 1] += temp_vec[1] * m[i];
            a[j * 3 + 2] += temp_vec[2] * m[i];
        }
    }
}

WIN32DLL_API void acceleration_massless(
    int objects_count,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    real G,
    real softening_length
)
{   
    real R_norm;
    real temp_value;
    real temp_vec[3];
    real R[3];
    
    // Empty the input array
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    int *restrict massive_indices = calloc(objects_count, sizeof(int));
    int *restrict massless_indices = calloc(objects_count, sizeof(int));
    int massive_objects_count = 0;
    int massless_objects_count = 0;
    for (int i = 0; i < objects_count; i++)
    {
        if (m[i] != 0.0)
        {
            massive_indices[massive_objects_count] = i;
            massive_objects_count++;
        }
        else
        {
            massless_indices[massless_objects_count] = i;
            massless_objects_count++;
        }
    }

    // Pairwise acceleration calculation for massive objects
    for (int i = 0; i < massive_objects_count; i++)
    {
        for (int j = i + 1; j < massive_objects_count; j++)
        {
            int idx_i = massive_indices[i];
            int idx_j = massive_indices[j];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2] + softening_length * softening_length);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[idx_i * 3 + 0] -= temp_vec[0] * m[j];
            a[idx_i * 3 + 1] -= temp_vec[1] * m[j];
            a[idx_i * 3 + 2] -= temp_vec[2] * m[j];
            a[idx_j * 3 + 0] += temp_vec[0] * m[i];
            a[idx_j * 3 + 1] += temp_vec[1] * m[i];
            a[idx_j * 3 + 2] += temp_vec[2] * m[i];
        }
    }

    // Acceleration calculation for massless objects
    for (int i = 0; i < massive_objects_count; i++)
    {
        for (int j = 0; j < massless_objects_count; j++)
        {
            int idx_i = massive_indices[i];
            int idx_j = massless_indices[j];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2] + softening_length * softening_length);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
            a[idx_j * 3 + 0] += temp_value * R[0] * m[i];
            a[idx_j * 3 + 1] += temp_value * R[1] * m[i];
            a[idx_j * 3 + 2] += temp_value * R[2] * m[i];
        }
    }

    free(massive_indices);
    free(massless_indices);
}


// // For debug
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
//     printf("Node (depth %d): is_leaf=%s, index=%d, mass=%.10f, center=(%.10f, %.10f, %.10f), width=%.2f\n",
//            depth, node->index >= 0 ? "true" : "false", node->index, node->total_mass,
//            node->center_of_mass[0], node->center_of_mass[1], node->center_of_mass[2],
//            node->box_width);

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
    // Find the width and center of the bounding box
    real center[3];
    real width;
    _calculate_bounding_box(objects_count, x, center, &width);

    // Construct the octree
    BarnesHutTreeNode *restrict root = malloc(sizeof(BarnesHutTreeNode));
    if (root != NULL) 
    {
        root->index = -1;
        root->total_mass = 0.0;
        root->center_of_mass[0] = center[0];
        root->center_of_mass[1] = center[1];
        root->center_of_mass[2] = center[2];
        root->box_width = width;
        root->children[0] = NULL;
        root->children[1] = NULL;
        root->children[2] = NULL;
        root->children[3] = NULL;
        root->children[4] = NULL;
        root->children[5] = NULL;
        root->children[6] = NULL;
        root->children[7] = NULL;
    }
    else
    {
        fprintf(stderr, "Error: Failed to allocate memory for the root node in acceleration_barnes_hut.\n");
        goto err_root_memory;
    }

    BarnesHutTreeNode *restrict child_node_pool = malloc(objects_count * sizeof(BarnesHutTreeNode));
    if (child_node_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for child nodes in acceleration_barnes_hut.\n");
        goto err_memory;
    }
    if (_barnes_hut_construct_octree(objects_count, x, m, width, child_node_pool, root) == 1)
    {
        goto err_memory;
    }

    // Calculate the center of mass
    if (_barnes_hut_compute_center_of_mass(root) == 1)
    {
        goto err_memory;
    }

    // Calculate the acceleration
    if (_barnes_hut_acceleration(barnes_hut_theta, objects_count, a, G, softening_length, root) == 1)
    {
        goto err_memory;
    }

    // Free the memory
    if (_barnes_hut_free_octree(child_node_pool, root) == 1)
    {
        fprintf(stderr, "Error: Failed to free memory of octree in acceleration_barnes_hut.\n");
        goto err_free_octree;
    }
    

    /* For debug / optimization
    clock_t start, end;

    start = clock();
    // Find the width of the bounding box
    real center[3];
    real width;
    _calculate_bounding_box(objects_count, x, center, &width);

    end = clock();
    printf("Time elapsed for finding the width: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    // Construct the octree
    BarnesHutTreeNode *restrict root = malloc(sizeof(BarnesHutTreeNode));
    if (root != NULL) 
    {
        root->index = -1;
        root->total_mass = 0.0;
        root->center_of_mass[0] = center[0];
        root->center_of_mass[1] = center[1];
        root->center_of_mass[2] = center[2];
        root->box_width = width;
        root->children[0] = NULL;
        root->children[1] = NULL;
        root->children[2] = NULL;
        root->children[3] = NULL;
        root->children[4] = NULL;
        root->children[5] = NULL;
        root->children[6] = NULL;
        root->children[7] = NULL;
    }
    else
    {
        fprintf(stderr, "Error: Failed to allocate memory for the root node in acceleration_barnes_hut.\n");
        goto err_root_memory;
    }

    BarnesHutTreeNode *restrict child_node_pool = malloc(objects_count * sizeof(BarnesHutTreeNode));
    if (child_node_pool == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for child nodes in acceleration_barnes_hut.\n");
        goto err_memory;
    }
    if (_barnes_hut_construct_octree(objects_count, x, m, width, child_node_pool, root) == 1)
    {
        goto err_memory;
    }
    // // Debug
    // visualize_octree_nodes(root, 0);
    // exit(1);

    end = clock();
    printf("Time elapsed for constructing the octree: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();

    // Calculate the center of mass
    if (_barnes_hut_compute_center_of_mass(root) == 1)
    {
        goto err_memory;
    }
    // // Debug
    // visualize_octree_nodes(root, 0);
    // exit(1);

    end = clock();
    printf("Time elapsed for computing the center of mass: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    // Calculate the acceleration
    if(_barnes_hut_acceleration(barnes_hut_theta, objects_count, a, G, softening_length, root) == 1)
    {
        goto err_memory;
    }

    end = clock();
    printf("Time elapsed for calculating the acceleration: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();

    // Free the memory
    if (_barnes_hut_free_octree(child_node_pool, root) == 1)
    {
        fprintf(stderr, "Error: Failed to free memory of octree in acceleration_barnes_hut.\n");
        goto err_free_octree;
    }
    end = clock();
    printf("Time elapsed for freeing the memory: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
    exit(0);
    */
    
    return;

err_memory:
    _barnes_hut_free_octree(child_node_pool, root);
err_free_octree:
err_root_memory:
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

    if (x > center_x) region |= 4;
    if (y > center_y) region |= 2;
    if (z > center_z) region |= 1;

    return region;
}

WIN32DLL_API int _barnes_hut_construct_octree(
    int objects_count,
    const real *restrict x,
    const real *restrict m,
    real width,
    BarnesHutTreeNode *restrict child_node_pool,
    BarnesHutTreeNode *restrict root
)
{
    for (int i = 0; i < objects_count; i++)
    {   
        int depth = 1;
        BarnesHutTreeNode *current_node = root;
        BarnesHutTreeNode *child_node = &child_node_pool[i];

        child_node->index = i;
        child_node->total_mass = m[i];
        child_node->center_of_mass[0] = x[i * 3 + 0];
        child_node->center_of_mass[1] = x[i * 3 + 1];
        child_node->center_of_mass[2] = x[i * 3 + 2];

        BarnesHutTreeNode *collider_leaf = NULL;
        
        while (true)
        {
            int region = _barnes_hut_check_region(
                x[i * 3 + 0],
                x[i * 3 + 1],
                x[i * 3 + 2],
                current_node->center_of_mass[0],
                current_node->center_of_mass[1],
                current_node->center_of_mass[2]
            );
            
            if (collider_leaf == NULL)
            {
                if (current_node->children[region] == NULL)
                {
                    current_node->children[region] = child_node;
                    break;
                }
                else if ((current_node->children[region]->index) >= 0)  // >= 0 means it is a leaf
                {
                    collider_leaf = current_node->children[region];
                    BarnesHutTreeNode *new_node = malloc(sizeof(BarnesHutTreeNode));
                    if (new_node == NULL)
                    {
                        goto err_memory;
                    }

                    new_node->index = -1;
                    new_node->total_mass = collider_leaf->total_mass + m[i];
                    new_node->box_width = width / (pow(2.0, depth));
                    new_node->children[0] = NULL;
                    new_node->children[1] = NULL;
                    new_node->children[2] = NULL;
                    new_node->children[3] = NULL;
                    new_node->children[4] = NULL;
                    new_node->children[5] = NULL;
                    new_node->children[6] = NULL;
                    new_node->children[7] = NULL;

                    // Calculate the center of the node when constructing the octree
                    // It will be replaced with the center of mass in the next step
                    if (region < 4)
                    {
                        new_node->center_of_mass[0] = current_node->center_of_mass[0] - width / (pow(2.0, depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[0] = current_node->center_of_mass[0] + width / (pow(2.0, depth + 1));
                    }

                    if (region == 0 || region == 1 || region == 4 || region == 5)
                    {
                        new_node->center_of_mass[1] = current_node->center_of_mass[1] - width / (pow(2.0, depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[1] = current_node->center_of_mass[1] + width / (pow(2.0, depth + 1));
                    }

                    if (region % 2 == 0)
                    {
                        new_node->center_of_mass[2] = current_node->center_of_mass[2] - width / (pow(2.0, depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[2] = current_node->center_of_mass[2] + width / (pow(2.0, depth + 1));
                    }

                    current_node->children[region] = new_node;
                    current_node = new_node;
                    depth++;
                }
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
                    current_node->children[region] = child_node;
                    current_node->children[collider_region] = collider_leaf;
                    break;
                }

                // Collide again
                else
                {
                    BarnesHutTreeNode *new_node = malloc(sizeof(BarnesHutTreeNode));
                    if (new_node == NULL)
                    {
                        goto err_memory;
                    }

                    new_node->index = -1;
                    new_node->total_mass = current_node->total_mass;
                    new_node->center_of_mass[0] = current_node->center_of_mass[0];
                    new_node->center_of_mass[1] = current_node->center_of_mass[1];
                    new_node->center_of_mass[2] = current_node->center_of_mass[2];
                    new_node->box_width = width / (pow(2.0, depth));
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
                    if (region < 4)
                    {
                        new_node->center_of_mass[0] = current_node->center_of_mass[0] - width / (pow(2.0, depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[0] = current_node->center_of_mass[0] + width / (pow(2.0, depth + 1));
                    }

                    if (region == 0 || region == 1 || region == 4 || region == 5)
                    {
                        new_node->center_of_mass[1] = current_node->center_of_mass[1] - width / (pow(2.0, depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[1] = current_node->center_of_mass[1] + width / (pow(2.0, depth + 1));
                    }
                    
                    if (region % 2 == 0)
                    {
                        new_node->center_of_mass[2] = current_node->center_of_mass[2] - width / (pow(2.0, depth + 1));
                    }
                    else
                    {
                        new_node->center_of_mass[2] = current_node->center_of_mass[2] + width / (pow(2.0, depth + 1));
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

WIN32DLL_API int _barnes_hut_compute_center_of_mass(BarnesHutTreeNode *restrict root)
{
    typedef struct BarnesHutStack
    {
        BarnesHutTreeNode *node;
        struct BarnesHutStack *last;
        real sum_of_mass_times_distance[3];
        int processed_region;
    } BarnesHutStack;
    
    BarnesHutTreeNode *current_node = root;
    BarnesHutStack *stack = malloc(sizeof(BarnesHutStack));
    if (stack == NULL)
    {
        goto err_memory;
    }
    stack->node = current_node;
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
            if (child_i == NULL)
            {
                stack->processed_region = i;
                continue;
            }
            else if (child_i->index >= 0)   // >= 0 means it is a leaf
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
                BarnesHutStack *new_item = malloc(sizeof(BarnesHutStack));
                if (new_item == NULL)
                {
                    goto err_memory;
                }
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
            real total_mass = current_node->total_mass;
            current_node->center_of_mass[0] = stack->sum_of_mass_times_distance[0] / total_mass;
            current_node->center_of_mass[1] = stack->sum_of_mass_times_distance[1] / total_mass;
            current_node->center_of_mass[2] = stack->sum_of_mass_times_distance[2] / total_mass;
        
            BarnesHutStack *parent = stack->last;
            free(stack);
            if (parent == NULL)
            {
                break;
            }

            BarnesHutTreeNode *parent_node = parent->node;
            parent->sum_of_mass_times_distance[0] += total_mass * current_node->center_of_mass[0];
            parent->sum_of_mass_times_distance[1] += total_mass * current_node->center_of_mass[1];
            parent->sum_of_mass_times_distance[2] += total_mass * current_node->center_of_mass[2];
            parent->processed_region++;

            stack = parent;
            current_node = parent_node;
        }
    }

    return 0;

err_memory:
    return 1;
}

WIN32DLL_API int _barnes_hut_acceleration(
    real theta,
    int objects_count,
    real *restrict a,
    real G,
    real softening_length,
    BarnesHutTreeNode *restrict root
)
{
    typedef struct BarnesHutStack
    {
        BarnesHutTreeNode *node;
        struct BarnesHutStack *last;
        int processed_region;
    } BarnesHutStack;

    BarnesHutTreeNode *current_acc_node = root;
    BarnesHutStack *acc_stack = malloc(sizeof(BarnesHutStack));
    if (acc_stack == NULL)
    {
        goto err_memory;
    }
    acc_stack->node = current_acc_node;
    acc_stack->last = NULL;
    acc_stack->processed_region = -1;

    // Empty the input array
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    BarnesHutTreeNode *current_acc_leaf = NULL;
    bool found_leaf;
    int acc_object_index = 0;
    BarnesHutTreeNode **acc_branch_nodes = calloc(objects_count, sizeof(BarnesHutTreeNode*));
    if (acc_branch_nodes == NULL)
    {
        goto err_memory;
    }
    int acc_branch_nodes_count;
    int i;
    while (true)
    {     
        // Find leaf
        found_leaf = false;
        acc_branch_nodes_count = 0;
        for (i = (acc_stack->processed_region + 1); i < 8; i++)
        {
            BarnesHutTreeNode *child_i = current_acc_node->children[i];
            if (child_i == NULL)
            {
                acc_stack->processed_region = i;
                continue;
            }
            else if (child_i->index >= 0)   // >= 0 means it is a leaf
            {
                current_acc_leaf = child_i;
                acc_object_index = current_acc_leaf->index;
                found_leaf = true;
                acc_stack->processed_region = i;
                break;
            }
            else
            {
                BarnesHutStack *new_item = malloc(sizeof(BarnesHutStack));
                if (new_item == NULL)
                {
                    goto err_memory;
                }
                new_item->node = child_i;
                new_item->last = acc_stack;
                new_item->processed_region = -1;

                acc_stack = new_item;
                current_acc_node = child_i;

                acc_branch_nodes[acc_branch_nodes_count] = current_acc_node;
                acc_branch_nodes_count++;
                break;
            }
        }
        
        if (found_leaf)
        {
            BarnesHutTreeNode *current_obj_node = root;
            BarnesHutStack *obj_stack = malloc(sizeof(BarnesHutStack));
            if (obj_stack == NULL)
            {
                goto err_memory;
            }
            obj_stack->node = root;
            obj_stack->last = NULL;
            obj_stack->processed_region = -1;

            // Calculate acceleration
            while (true)
            {
                for (int j = (obj_stack->processed_region + 1); j < 8; j++)
                {
                    BarnesHutTreeNode *child_j = current_obj_node->children[j];
                    if (child_j == NULL)
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

                        bool is_acc_branch_node = false;
                        for (int k = 0; k < acc_branch_nodes_count; k++)
                        {
                            if (acc_branch_nodes[k] == child_j)
                            {
                                is_acc_branch_node = true;
                                break;
                            }
                        }

                        if (
                            ((child_j->box_width / R_norm) < theta && (!is_acc_branch_node)) ||
                            (child_j->index >= 0)   // >= 0 means it is a leaf
                        )
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
                            BarnesHutStack *new_item = malloc(sizeof(BarnesHutStack));
                            if (new_item == NULL)
                            {
                                goto err_memory;
                            }
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
                    BarnesHutStack *parent = obj_stack->last;
                    free(obj_stack);
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
            BarnesHutStack *parent = acc_stack->last;
            free(acc_stack);
            if (parent == NULL)
            {
                break;
            }
            parent->processed_region++;
            acc_stack = parent;
            current_acc_node = parent->node;
        }
    }

    free(acc_branch_nodes);
    return 0;

err_memory:
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

WIN32DLL_API int _barnes_hut_free_octree(
    BarnesHutTreeNode *restrict child_node_pool,
    BarnesHutTreeNode *restrict root
)
{
    typedef struct BarnesHutStack
    {
        BarnesHutTreeNode *node;
        struct BarnesHutStack *last;
        int processed_region;
    } BarnesHutStack;
    
    BarnesHutTreeNode *current_node = root;
    BarnesHutStack *stack = malloc(sizeof(BarnesHutStack));
    if (stack == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate memory for the stack in _barnes_hut_free_octree.\n");
        goto err_memory;
    }
    stack->node = current_node;
    stack->last = NULL;
    stack->processed_region = -1;

    while (true)
    {   
        for (int i = (stack->processed_region + 1); i < 8; i++)
        {
            BarnesHutTreeNode *child_i = current_node->children[i];
            if (child_i == NULL || child_i->index >= 0) // >= 0 means it is a leaf
            {
                stack->processed_region = i;
                continue;
            }
            else
            {
                BarnesHutStack *new_item = malloc(sizeof(BarnesHutStack));
                if (new_item == NULL)
                {
                    goto err_memory;
                }
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
            BarnesHutStack *parent = stack->last;
            free(current_node);
            free(stack);
            if (parent == NULL)
            {
                break;
            }

            parent->processed_region++;
            current_node = parent->node;
            stack = parent;
        }
    }
    free(child_node_pool);

    return 0;

err_memory:
    return 1;
}
