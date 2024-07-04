#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

void ias15(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real (*restrict a)[3], 
    const real *restrict m, 
    real G,    
    int dim_nodes,
    const real *restrict nodes,
    const real *restrict aux_c, 
    const real *restrict aux_r,
    real *restrict aux_b0,
    real *restrict aux_b,
    real *restrict aux_g,
    real *restrict aux_e, 
    real *restrict t, 
    real *restrict dt, 
    real expected_time_scale, 
    real tolerance,
    real tolerance_pc,
    real safety_fac,
    real exponent,
    int *restrict ias15_refine_flag,
    int max_iteration,
    int min_iteration
);
void ias15_step(
    int objects_count,
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    const real *restrict m,
    real G,
    real *restrict t,
    real *restrict dt,
    real expected_time_scale,
    int dim_nodes,
    int dim_nodes_minus_1,
    int dim_nodes_minus_2,
    const real *restrict nodes,
    real *restrict aux_b0,
    real *restrict aux_b,
    const real *restrict aux_c,
    real *restrict aux_e,
    real *restrict aux_g,
    const real *restrict aux_r,
    real tolerance,
    real tolerance_pc,
    real exponent,
    real safety_fac,
    int *restrict ias15_refine_flag,
    real *restrict aux_a,
    real (*restrict x)[3],
    real (*restrict v)[3],
    real (*restrict a)[3],
    real *restrict delta_b7,
    real *restrict F,
    real *restrict delta_aux_b
);
void ias15_approx_pos(
    int objects_count,
    real (*restrict x)[3],
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real node,
    real *restrict aux_b,
    real dt
);
void ias15_approx_vel(
    int objects_count,
    real (*restrict v)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real node,
    real *restrict aux_b,
    real dt
);
void ias15_compute_aux_b(
    int objects_count,
    int dim_nodes_minus_1,
    real *restrict aux_b,
    const real *restrict aux_g,
    const real *restrict aux_c,
    int i
);
void ias15_compute_aux_g(
    int objects_count,
    int dim_nodes,
    real *restrict aux_g,
    const real *restrict aux_r,
    const real *restrict aux_a,
    int i,
    real *restrict F
);
void ias15_refine_aux_b(
    int objects_count,
    int dim_nodes_minus_1,
    real *restrict aux_b,
    real *restrict aux_e,
    real *restrict delta_aux_b,
    real dt,
    real dt_new,
    int ias15_refine_flag
);


WIN32DLL_API void ias15(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real (*restrict a)[3], 
    const real *restrict m, 
    real G,    
    int dim_nodes,
    const real *restrict nodes,
    const real *restrict aux_c, 
    const real *restrict aux_r,
    real *restrict aux_b0,
    real *restrict aux_b,
    real *restrict aux_g,
    real *restrict aux_e, 
    real *restrict t, 
    real *restrict dt, 
    real expected_time_scale, 
    real tolerance,
    real tolerance_pc,
    real safety_fac,
    real exponent,
    int *restrict ias15_refine_flag,
    int max_iteration,
    int min_iteration
)
{
    real t0 = *t;

    int dim_nodes_minus_1 = dim_nodes - 1;
    int dim_nodes_minus_2 = dim_nodes - 2;

    // Arrays for ias15_step
    real *aux_a = malloc(dim_nodes * objects_count * 3 * sizeof(real));
    real (*x_step)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*v_step)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a_step)[3] = malloc(objects_count * 3 * sizeof(real));
    real *delta_b7 = malloc(objects_count * 3 * sizeof(real));

    // Arrays for compute aux_g
    real *F = malloc(8 * objects_count * 3 * sizeof(real));

    // Array for refine aux_b
    real *delta_aux_b = malloc(dim_nodes_minus_1 * objects_count * 3 * sizeof(real));

    for (int i = 0; i < max_iteration; i++)
    {
        ias15_step(
            objects_count,
            x,
            v,
            a,
            m,
            G,
            t,
            dt,
            expected_time_scale,
            dim_nodes,
            dim_nodes_minus_1,
            dim_nodes_minus_2,
            nodes,
            aux_b0,
            aux_b,
            aux_c,
            aux_e,
            aux_g,
            aux_r,
            tolerance,
            tolerance_pc,
            exponent,
            safety_fac,
            ias15_refine_flag,
            aux_a,
            x_step,
            v_step,
            a_step,
            delta_b7,
            F,
            delta_aux_b
        );

        if (i >= min_iteration && *t > (t0 + expected_time_scale * 1e-5))
        {
            free(aux_a);
            free(x_step);
            free(v_step);
            free(a_step);
            free(delta_b7); 
            free(F);
            free(delta_aux_b);
            break;
        }
    }
}

// Advance IAS15 for one step
WIN32DLL_API void ias15_step(
    int objects_count,
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    const real *restrict m,
    real G,
    real *restrict t,
    real *restrict dt,
    real expected_time_scale,
    int dim_nodes,
    int dim_nodes_minus_1,
    int dim_nodes_minus_2,
    const real *restrict nodes,
    real *restrict aux_b0,
    real *restrict aux_b,
    const real *restrict aux_c,
    real *restrict aux_e,
    real *restrict aux_g,
    const real *restrict aux_r,
    real tolerance,
    real tolerance_pc,
    real exponent,
    real safety_fac,
    int *restrict ias15_refine_flag,
    real *restrict aux_a,
    real (*restrict x)[3],
    real (*restrict v)[3],
    real (*restrict a)[3],
    real *restrict delta_b7,
    real *restrict F,
    real *restrict delta_aux_b
)
{
    real error, error_b7, dt_new;
    // Main Loop
    int ias15_integrate_flag = 0; 
    while (1)
    {   
        // Loop for predictor-corrector algorithm
        // 12 = max iterations
        for (int temp = 0; temp < 12; temp++)
        {
            for (int i = 0; i < dim_nodes; i++)
            {
                // Estimate position and velocity with current aux_b and nodes
                ias15_approx_pos(objects_count, x, x0, v0, a0, nodes[i], aux_b, *dt);
                ias15_approx_vel(objects_count, v, v0, a0, nodes[i], aux_b, *dt);

                // Evaluate force function and store result
                acceleration(objects_count, x, a, m, G);
                memcpy(&aux_a[i * objects_count * 3], a, objects_count * 3 * sizeof(real));
                
                ias15_compute_aux_g(objects_count, dim_nodes, aux_g, aux_r, aux_a, i, F);
                ias15_compute_aux_b(objects_count, dim_nodes_minus_1, aux_b, aux_g, aux_c, i);
            }

            // Estimate convergence
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    delta_b7[i * 3 + j] = aux_b[dim_nodes_minus_2 * objects_count * 3 + i * 3 + j] - aux_b0[dim_nodes_minus_2 * objects_count * 3 + i * 3 + j];
                }
            }
            memcpy(aux_b0, aux_b, dim_nodes_minus_1 * objects_count * 3 * sizeof(real));
            if ((abs_max_vec(delta_b7, objects_count * 3) / abs_max_vec(&aux_a[dim_nodes_minus_1 * objects_count * 3], objects_count * 3)) < tolerance_pc)
            {
                break;
            }
        }

        // Advance step
        ias15_approx_pos(objects_count, x, x0, v0, a0, 1.0, aux_b, *dt);
        ias15_approx_vel(objects_count, v, v0, a0, 1.0, aux_b, *dt);
        acceleration(objects_count, x, a, m, G);

        // Estimate relative error
        error_b7 = abs_max_vec(&aux_b[dim_nodes_minus_2 * objects_count * 3], objects_count * 3) / abs_max_vec_array(a, objects_count);
        error = pow((error_b7 / tolerance), exponent);
        
        // Step-size for the next step
        if (error != 0)
        {
            dt_new = *dt / error;
        }
        else
        {
            dt_new = *dt;
        }

        // Accept the step
        if (error <= 1 || *dt == expected_time_scale * 1e-12)
        {
            // Report accepted step
            ias15_integrate_flag = 1;
            *t += *dt;

            ias15_refine_aux_b(objects_count, dim_nodes_minus_1, aux_b, aux_e, delta_aux_b, *dt, dt_new, *ias15_refine_flag);
            *ias15_refine_flag = 1;
        }

        // Step size for the next iteration
        if (dt_new > (*dt / safety_fac))
        {
            *dt = *dt / safety_fac;
        }
        else if (dt_new < *dt * safety_fac)
        {
            *dt = *dt * safety_fac;
        }
        else
        {
            *dt = dt_new;
        }

        if (dt_new / expected_time_scale < 1e-12)
        {
            *dt = expected_time_scale * 1e-12;
        }

        // Exit 
        if (ias15_integrate_flag > 0)
        {   
            memcpy(x0, x, objects_count * 3 * sizeof(real));
            memcpy(v0, v, objects_count * 3 * sizeof(real));
            memcpy(a0, a, objects_count * 3 * sizeof(real));      
            break;    
        }
    }
}

WIN32DLL_API void ias15_approx_pos(
    int objects_count,
    real (*restrict x)[3],
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real node,
    real *restrict aux_b,
    real dt
)
{   
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            /*
            *   Warning: Combining both statements would increase floating point error
            *            e.g. x[j][k] = x0[j][k] + ...   (WRONG)
            */

            x[j][k] = x0[j][k];
            x[j][k] += dt * node * (
                v0[j][k]
                + dt
                * node
                * (
                    a0[j][k]
                    + node 
                    * (
                        aux_b[0 * objects_count * 3 + j * 3 + k] / 3.0
                        + node
                        * (
                            aux_b[1 * objects_count * 3 + j * 3 + k] / 6.0
                            + node
                            * (
                                aux_b[2 * objects_count * 3 + j * 3 + k] / 10.0
                                + node
                                * (
                                    aux_b[3 * objects_count * 3 + j * 3 + k] / 15.0
                                    + node
                                    * (
                                        aux_b[4 * objects_count * 3 + j * 3 + k] / 21.0
                                        + node * (
                                            aux_b[5 * objects_count * 3 + j * 3 + k] / 28.0 
                                            + node * aux_b[6 * objects_count * 3 + j * 3 + k] / 36.0
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
                / 2.0
            );
        }
    }
}

WIN32DLL_API void ias15_approx_vel(
    int objects_count,
    real (*restrict v)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real node,
    real *restrict aux_b,
    real dt
)
{
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            /*
            *   Warning: Combining both statements would increase floating point error
            *            e.g. v[j][k] = v0[j][k] + ...   (WRONG)
            */

            v[j][k] = v0[j][k];
            v[j][k] += dt * node * (
                a0[j][k]
                + node
                * (
                    aux_b[0 * objects_count * 3 + j * 3 + k] / 2.0
                    + node
                    * (
                        aux_b[1 * objects_count * 3 + j * 3 + k] / 3.0
                        + node
                        * (
                            aux_b[2 * objects_count * 3 + j * 3 + k] / 4.0
                            + node
                            * (
                                aux_b[3 * objects_count * 3 + j * 3 + k] / 5.0
                                + node
                                * (
                                    aux_b[4 * objects_count * 3 + j * 3 + k] / 6.0
                                    + node * (
                                        aux_b[5 * objects_count * 3 + j * 3 + k] / 7.0 
                                        + node * aux_b[6 * objects_count * 3 + j * 3 + k] / 8.0
                                    )
                                )
                            )
                        )
                    )
                )
            );
        }
    }
}

// Calculate the auxiliary coefficients b for IAS15
WIN32DLL_API void ias15_compute_aux_b(
    int objects_count,
    int dim_nodes_minus_1,
    real *restrict aux_b,
    const real *restrict aux_g,
    const real *restrict aux_c,
    int i
)
{
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            if (i >= 1) {
                aux_b[0 * objects_count * 3 + j * 3 + k] = (
                    aux_c[0 * dim_nodes_minus_1 + 0] * aux_g[0 * objects_count * 3 + j * 3 + k]
                    + aux_c[1 * dim_nodes_minus_1 + 0] * aux_g[1 * objects_count * 3 + j * 3 + k]
                    + aux_c[2 * dim_nodes_minus_1 + 0] * aux_g[2 * objects_count * 3 + j * 3 + k]
                    + aux_c[3 * dim_nodes_minus_1 + 0] * aux_g[3 * objects_count * 3 + j * 3 + k]
                    + aux_c[4 * dim_nodes_minus_1 + 0] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 0] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 0] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 2) {
                aux_b[1 * objects_count * 3 + j * 3 + k] = (
                    aux_c[1 * dim_nodes_minus_1 + 1] * aux_g[1 * objects_count * 3 + j * 3 + k]
                    + aux_c[2 * dim_nodes_minus_1 + 1] * aux_g[2 * objects_count * 3 + j * 3 + k]
                    + aux_c[3 * dim_nodes_minus_1 + 1] * aux_g[3 * objects_count * 3 + j * 3 + k]
                    + aux_c[4 * dim_nodes_minus_1 + 1] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 1] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 1] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 3) {
                aux_b[2 * objects_count * 3 + j * 3 + k] = (
                    aux_c[2 * dim_nodes_minus_1 + 2] * aux_g[2 * objects_count * 3 + j * 3 + k]
                    + aux_c[3 * dim_nodes_minus_1 + 2] * aux_g[3 * objects_count * 3 + j * 3 + k]
                    + aux_c[4 * dim_nodes_minus_1 + 2] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 2] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 2] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 4) {
                aux_b[3 * objects_count * 3 + j * 3 + k] = (
                    aux_c[3 * dim_nodes_minus_1 + 3] * aux_g[3 * objects_count * 3 + j * 3 + k]
                    + aux_c[4 * dim_nodes_minus_1 + 3] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 3] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 3] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 5)
            {
                aux_b[4 * objects_count * 3 + j * 3 + k] = (
                    aux_c[4 * dim_nodes_minus_1 + 4] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 4] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 4] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 6)
            {
                aux_b[5 * objects_count * 3 + j * 3 + k] = (
                    aux_c[5 * dim_nodes_minus_1 + 5] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 5] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 7)
            {
                aux_b[6 * objects_count * 3 + j * 3 + k] = (
                    aux_c[6 * dim_nodes_minus_1 + 6] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }
        }
    }
}

WIN32DLL_API void ias15_compute_aux_g(
    int objects_count,
    int dim_nodes,
    real *restrict aux_g,
    const real *restrict aux_r,
    const real *restrict aux_a,
    int i,
    real *restrict F
)
{
    // Retrieve required accelerations
    // F is allocated in IAS15
    for (int j = 0; j <= i; j++)
    {
        memcpy(&F[j * objects_count * 3], &aux_a[j * objects_count * 3], objects_count * 3 * sizeof(real));
    }
    
    // Update aux_g
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            if (i >= 1)
            {
                aux_g[0 * objects_count * 3 + j * 3 + k] = (
                    (F[1 * objects_count * 3 + j * 3 + k] - F[0 * objects_count * 3 + j * 3 + k]) * aux_r[1 * dim_nodes + 0]
                );
            }
            else
            {
                continue;
            }
                
            if (i >= 2)
            {
                aux_g[1 * objects_count * 3 + j * 3 + k] = (
                    ((F[2 * objects_count * 3 + j * 3 + k] - F[0 * objects_count * 3 + j * 3 + k]) * aux_r[2 * dim_nodes + 0] 
                    - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[2 * dim_nodes + 1]
                );
            }
            else 
            {
                continue;
            }

            if (i >= 3)
            {
                aux_g[2 * objects_count * 3 + j * 3 + k] = (
                    ((F[3 * objects_count * 3 + j * 3 + k] - F[0 * objects_count * 3 + j * 3 + k]) * aux_r[3 * dim_nodes + 0] 
                    - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[3 * dim_nodes + 1] 
                    - aux_g[1 * objects_count * 3 + j * 3 + k]
                ) * aux_r[3 * dim_nodes + 2];
            }
            else
            {
                continue;
            }
                
            if (i >= 4)
            {
                aux_g[3 * objects_count * 3 + j * 3 + k] = (
                    (((F[4 * objects_count * 3 + j * 3 + k] - F[0 * objects_count * 3 + j * 3 + k]) * aux_r[4 * dim_nodes + 0] 
                    - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[4 * dim_nodes + 1] - aux_g[1 * objects_count * 3 + j * 3 + k])
                    * aux_r[4 * dim_nodes + 2]
                    - aux_g[2 * objects_count * 3 + j * 3 + k]
                ) * aux_r[4 * dim_nodes + 3];
            }
            else
            {
                continue;
            }

            if (i >= 5)
            {
                aux_g[4 * objects_count * 3 + j * 3 + k] = (
                    (
                        (((F[5 * objects_count * 3 + j * 3 + k] - F[0 * objects_count * 3 + j * 3 + k]) * aux_r[5 * dim_nodes + 0] 
                        - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[5 * dim_nodes + 1] 
                        - aux_g[1 * objects_count * 3 + j * 3 + k])
                        * aux_r[5 * dim_nodes + 2]
                        - aux_g[2 * objects_count * 3 + j * 3 + k]
                    )
                    * aux_r[5 * dim_nodes + 3]
                    - aux_g[3 * objects_count * 3 + j * 3 + k]
                ) * aux_r[5 * dim_nodes + 4];
            }
            else
            {
                continue;
            }

            if (i >= 6)
            {
                aux_g[5 * objects_count * 3 + j * 3 + k] = (
                    (
                        (
                            (((F[6 * objects_count * 3 + j * 3 + k] - F[0 * objects_count * 3 + j * 3 + k]) * aux_r[6 * dim_nodes + 0] 
                            - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[6 * dim_nodes + 1] 
                            - aux_g[1 * objects_count * 3 + j * 3 + k])
                            * aux_r[6 * dim_nodes + 2]
                            - aux_g[2 * objects_count * 3 + j * 3 + k]
                        )
                        * aux_r[6 * dim_nodes + 3]
                        - aux_g[3 * objects_count * 3 + j * 3 + k]
                    )
                    * aux_r[6 * dim_nodes + 4]
                    - aux_g[4 * objects_count * 3 + j * 3 + k]
                ) * aux_r[6 * dim_nodes + 5];
            }
            else
            {
                continue;
            }

            if (i >= 7)
            {
                aux_g[6 * objects_count * 3 + j * 3 + k] = (
                    (
                        (
                            (
                                (((F[7 * objects_count * 3 + j * 3 + k] - F[0 * objects_count * 3 + j * 3 + k]) * aux_r[7 * dim_nodes + 0] - aux_g[0 * objects_count * 3 + j * 3 + k]) 
                                * aux_r[7 * dim_nodes + 1] 
                                - aux_g[1 * objects_count * 3 + j * 3 + k])
                                * aux_r[7 * dim_nodes + 2]
                                - aux_g[2 * objects_count * 3 + j * 3 + k]
                            )
                            * aux_r[7 * dim_nodes + 3]
                            - aux_g[3 * objects_count * 3 + j * 3 + k]
                        )
                        * aux_r[7 * dim_nodes + 4]
                        - aux_g[4 * objects_count * 3 + j * 3 + k]
                    )
                    * aux_r[7 * dim_nodes + 5]
                    - aux_g[5 * objects_count * 3 + j * 3 + k]
                ) * aux_r[7 * dim_nodes + 6];                
            }
            else
            {
                continue;
            }
        }
    }
}

WIN32DLL_API void ias15_refine_aux_b(
    int objects_count,
    int dim_nodes_minus_1,
    real *restrict aux_b,
    real *restrict aux_e,
    real *restrict delta_aux_b,
    real dt,
    real dt_new,
    int ias15_refine_flag
)
{
    if (ias15_refine_flag != 0)
    {
        for (int i = 0; i < dim_nodes_minus_1; i++)
        {
            for (int j = 0; j < objects_count; j++)
            {
                delta_aux_b[i * objects_count * 3 + j * 3 + 0] = (
                    aux_b[i * objects_count * 3 + j * 3 + 0] - aux_e[i * objects_count * 3 + j * 3 + 0]
                );
                delta_aux_b[i * objects_count * 3 + j * 3 + 1] = (
                    aux_b[i * objects_count * 3 + j * 3 + 1] - aux_e[i * objects_count * 3 + j * 3 + 1]
                );
                delta_aux_b[i * objects_count * 3 + j * 3 + 2] = (
                    aux_b[i * objects_count * 3 + j * 3 + 2] - aux_e[i * objects_count * 3 + j * 3 + 2]
                );
            }
        }
    }
    else
    {
        // Empty delta_aux_b
        for (int i = 0; i < dim_nodes_minus_1 * objects_count; i++)
        {
            delta_aux_b[i * 3 + 0] = 0.0;
            delta_aux_b[i * 3 + 1] = 0.0;
            delta_aux_b[i * 3 + 2] = 0.0;
        }
    }

    real q = dt_new / dt;
    real q2 = q * q, q3 = q2 * q, q4 = q3 * q, q5 = q4 * q, q6 = q5 * q, q7 = q6 * q;
    
    for (int j = 0; j < objects_count; j++) 
    {
        for (int k = 0; k < 3; k++) 
        {
            aux_e[0 * objects_count * 3 + j * 3 + k] = q * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 7.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 6.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] * 5.0
                + aux_b[3 * objects_count * 3 + j * 3 + k] * 4.0
                + aux_b[2 * objects_count * 3 + j * 3 + k] * 3.0
                + aux_b[1 * objects_count * 3 + j * 3 + k] * 2.0
                + aux_b[0 * objects_count * 3 + j * 3 + k]
            );

            aux_e[1 * objects_count * 3 + j * 3 + k] = q2 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 21.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 15.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] * 10.0
                + aux_b[3 * objects_count * 3 + j * 3 + k] * 6.0
                + aux_b[2 * objects_count * 3 + j * 3 + k] * 3.0
                + aux_b[1 * objects_count * 3 + j * 3 + k]
            );

            aux_e[2 * objects_count * 3 + j * 3 + k] = q3 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 35.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 20.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] * 10.0
                + aux_b[3 * objects_count * 3 + j * 3 + k] * 4.0
                + aux_b[2 * objects_count * 3 + j * 3 + k]
            );

            aux_e[3 * objects_count * 3 + j * 3 + k] = q4 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 35.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 15.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] * 5.0
                + aux_b[3 * objects_count * 3 + j * 3 + k]
            );

            aux_e[4 * objects_count * 3 + j * 3 + k] = q5 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 21.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 6.0
                + aux_b[4 * objects_count * 3 + j * 3 + k]
            );

            aux_e[5 * objects_count * 3 + j * 3 + k] = q6 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 7.0
                + aux_b[5 * objects_count * 3 + j * 3 + k]
            );

            aux_e[6 * objects_count * 3 + j * 3 + k] = q7 * aux_b[6 * objects_count * 3 + j * 3 + k];
        }
    }

    for (int i = 0; i < dim_nodes_minus_1; i++)
    {
        for (int j = 0; j < objects_count; j++)
        {
            aux_b[i * objects_count * 3 + j * 3 + 0] = (
                    aux_e[i * objects_count * 3 + j * 3 + 0] + delta_aux_b[i * objects_count * 3 + j * 3 + 0]
                );
            aux_b[i * objects_count * 3 + j * 3 + 1] = (
                aux_e[i * objects_count * 3 + j * 3 + 1] + delta_aux_b[i * objects_count * 3 + j * 3 + 1]
            );
            aux_b[i * objects_count * 3 + j * 3 + 2] = (
                aux_e[i * objects_count * 3 + j * 3 + 2] + delta_aux_b[i * objects_count * 3 + j * 3 + 2]
            );
        }
    }
}
