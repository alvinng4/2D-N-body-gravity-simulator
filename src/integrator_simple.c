/**
 * \file integrator_simple.c
 * \author Ching Yin Ng
 * \brief Function definitions for simple integrators
 * 
 * This file contains the function definitions for simple integrators,
 * including Euler, Euler-Cromer, Runge-Kutta 4th order (RK4),
 * and Leapfrog.
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acceleration.h"
#include "error.h"
#include "gravity_sim.h"
#include "storing.h"

WIN32DLL_API int euler(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
)
{
    /* Declare variables */
    int return_code;

    real *restrict x = system->x;
    real *restrict v = system->v;
    const int objects_count = system->objects_count;

    const real dt = integrator_param->dt;
    const int64 n_steps = simulation_param->n_steps_;

    const int storing_freq = storing_param->storing_freq;

    /* Allocate memory */
    real *restrict x_0 = malloc(objects_count * 3 * sizeof(real));
    real *restrict v_0 = malloc(objects_count * 3 * sizeof(real));
    real *restrict a = malloc(objects_count * 3 * sizeof(real));

    // Compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    // Check if memory allocation is successful
    if (!x_0 || !v_0 || !a || !x_err_comp_sum || !v_err_comp_sum)
    {
        return_code = ERROR_EULER_MEMORY_ALLOC;
        goto err_memory;
    }

    /* Main Loop */
    for (int64 count = 1; count <= n_steps; count++)
    {   
        memcpy(x_0, x, objects_count * 3 * sizeof(real));
        memcpy(v_0, v, objects_count * 3 * sizeof(real));

        /* Compute acceleration */
        return_code = acceleration(
            a,
            system,
            acceleration_param
        );
        if (return_code != SUCCESS)
        {
            goto err_acceleration;
        }

        /* Update step */
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++) 
            {
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;

                x[i * 3 + j] = x_0[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                v[i * 3 + j] = v_0[i * 3 + j] + v_err_comp_sum[i * 3 + j];

                x_err_comp_sum[i * 3 + j] += x_0[i * 3 + j] - x[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += v_0[i * 3 + j] - v[i * 3 + j];
            }
        }
        *(simulation_status->t) = count * dt;

        /* Store solution */
        if (count % storing_freq == 0)
        {
            return_code = store_solution_step(
                storing_param,
                system,
                simulation_status,
                solutions
            );
            if (return_code != SUCCESS)
            {
                goto err_store_solution;
            }
        }

        /* Check user interrupt */
        if (*(settings->is_exit))
        {
            return_code = ERROR_USER_INTERRUPT;
            goto err_user_interrupt;
        }
    }

    return SUCCESS;

err_user_interrupt:
err_store_solution:
err_acceleration:
err_memory:
    free(x_0);
    free(v_0);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    return return_code;
}

WIN32DLL_API int euler_cromer(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
)
{
    /* Declare variables */
    int return_code;

    real *restrict x = system->x;
    real *restrict v = system->v;
    const int objects_count = system->objects_count;

    const real dt = integrator_param->dt;
    const int64 n_steps = simulation_param->n_steps_;

    const int storing_freq = storing_param->storing_freq;

    /* Allocate memory */
    real *restrict x_0 = malloc(objects_count * 3 * sizeof(real));
    real *restrict v_0 = malloc(objects_count * 3 * sizeof(real));
    real *restrict a = malloc(objects_count * 3 * sizeof(real));
    
    // Compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    // Check if memory allocation is successful
    if (!x_0 || !v_0 || !a || !x_err_comp_sum || !v_err_comp_sum)
    {
        return_code = ERROR_EULER_CROMER_MEMORY_ALLOC;
        goto err_memory;
    }

    /* Main Loop */
    for (int64 count = 1; count <= n_steps; count++)
    {   
        memcpy(x_0, x, objects_count * 3 * sizeof(real));
        memcpy(v_0, v, objects_count * 3 * sizeof(real));

        /* Compute acceleration */
        return_code = acceleration(
            a,
            system,
            acceleration_param
        );
        if (return_code != SUCCESS)
        {
            goto err_acceleration;
        }

        /* Update step */
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++) 
            {
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;
                v[i * 3 + j] = v_0[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += v_0[i * 3 + j] - v[i * 3 + j];

                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                x[i * 3 + j] = x_0[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += x_0[i * 3 + j] - x[i * 3 + j];
            }
        }
        *(simulation_status->t) = count * dt;

        /* Store solution */
        if (count % storing_freq == 0)
        {
            return_code = store_solution_step(
                storing_param,
                system,
                simulation_status,
                solutions
            );
            if (return_code != SUCCESS)
            {
                goto err_store_solution;
            }
        }

        /* Check user interrupt */
        if (*(settings->is_exit))
        {
            return_code = ERROR_USER_INTERRUPT;
            goto err_user_interrupt;
        }
    }

    return SUCCESS;

err_user_interrupt:
err_store_solution:
err_acceleration:
err_memory:
    free(x_0);
    free(v_0);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    return return_code;
}

WIN32DLL_API int rk4(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
)
{
    /* Declare variables */
    int return_code;

    real *restrict x = system->x;
    real *restrict v = system->v;
    const int objects_count = system->objects_count;

    const real dt = integrator_param->dt;
    const int64 n_steps = simulation_param->n_steps_;

    const int storing_freq = storing_param->storing_freq;

    /* Allocate memory */
    real *restrict x_0 = malloc(objects_count * 3 * sizeof(real));
    real *restrict v_0 = malloc(objects_count * 3 * sizeof(real));
    real *restrict vk1 = malloc(objects_count * 3 * sizeof(real));
    real *restrict vk2 = malloc(objects_count * 3 * sizeof(real));
    real *restrict vk3 = malloc(objects_count * 3 * sizeof(real));
    real *restrict vk4 = malloc(objects_count * 3 * sizeof(real));
    real *restrict xk1 = malloc(objects_count * 3 * sizeof(real));
    real *restrict xk2 = malloc(objects_count * 3 * sizeof(real));
    real *restrict xk3 = malloc(objects_count * 3 * sizeof(real));
    real *restrict xk4 = malloc(objects_count * 3 * sizeof(real));
    
    // Compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    // Check if memory allocation is successful
    if (!x_0 
        || !v_0 
        || !vk1 
        || !vk2 
        || !vk3 
        || !vk4 
        || !xk1 
        || !xk2 
        || !xk3 
        || !xk4 
        || !x_err_comp_sum 
        || !v_err_comp_sum
    )
    {
        return_code = ERROR_RK4_MEMORY_ALLOC;
        goto err_memory;
    }

    /* Main Loop */
    for (int64 count = 1; count <= n_steps; count++)
    {
        // Store current state
        memcpy(x_0, x, objects_count * 3 * sizeof(real));
        memcpy(v_0, v, objects_count * 3 * sizeof(real));

        /* Compute xk1 and vk1 */
        return_code = acceleration(
            vk1,
            system,
            acceleration_param
        );
        if (return_code != SUCCESS)
        {
            goto err_acceleration;
        }
        memcpy(xk1, v, objects_count * 3 * sizeof(real));

        /* Compute xk2 and vk2 */
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = x_0[i * 3 + j] + 0.5 * xk1[i * 3 + j] * dt;
                v[i * 3 + j] = v_0[i * 3 + j] + 0.5 * vk1[i * 3 + j] * dt;
            }
        }
        return_code = acceleration(
            vk2,
            system,
            acceleration_param
        );
        if (return_code != SUCCESS)
        {
            goto err_acceleration;
        }
        memcpy(xk2, v, objects_count * 3 * sizeof(real));

        /* Compute xk3 and vk3 */
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = x_0[i * 3 + j] + 0.5 * xk2[i * 3 + j] * dt;
                v[i * 3 + j] = v_0[i * 3 + j] + 0.5 * vk2[i * 3 + j] * dt;
            }
        }
        return_code = acceleration(
            vk3,
            system,
            acceleration_param
        );
        memcpy(xk3, v, objects_count * 3 * sizeof(real));

        /* Compute xk4 and vk4 */
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = x_0[i * 3 + j] + xk3[i * 3 + j] * dt;
                v[i * 3 + j] = v_0[i * 3 + j] + vk3[i * 3 + j] * dt;
            }
        }
        return_code = acceleration(
            vk4,
            system,
            acceleration_param
        );
        memcpy(xk4, v, objects_count * 3 * sizeof(real));

        /* Update step */
        memcpy(v, v_0, objects_count * 3 * sizeof(real));
        memcpy(x, x_0, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                v_err_comp_sum[i * 3 + j] += (vk1[i * 3 + j] + 2 * vk2[i * 3 + j] + 2 * vk3[i * 3 + j] + vk4[i * 3 + j]) * dt / 6.0;
                x_err_comp_sum[i * 3 + j] += (xk1[i * 3 + j] + 2 * xk2[i * 3 + j] + 2 * xk3[i * 3 + j] + xk4[i * 3 + j]) * dt / 6.0;

                v[i * 3 + j] = v_0[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                x[i * 3 + j] = x_0[i * 3 + j] + x_err_comp_sum[i * 3 + j];

                v_err_comp_sum[i * 3 + j] += v_0[i * 3 + j] - v[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += x_0[i * 3 + j] - x[i * 3 + j];
            }
        }
        *(simulation_status->t) = count * dt;

        /* Store solution */
        if (count % storing_freq == 0)
        {
            return_code = store_solution_step(
                storing_param,
                system,
                simulation_status,
                solutions
            );
            if (return_code != SUCCESS)
            {
                goto err_store_solution;
            }
        }

        /* Check user interrupt */
        if (*(settings->is_exit))
        {
            return_code = ERROR_USER_INTERRUPT;
            goto err_user_interrupt;
        }
    }

    return SUCCESS;

err_user_interrupt:
err_store_solution:
err_acceleration:
err_memory:
    free(x_0);
    free(v_0);
    free(vk1);
    free(vk2);
    free(vk3);
    free(vk4);
    free(xk1);
    free(xk2);
    free(xk3);
    free(xk4);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    return return_code;
}

WIN32DLL_API int leapfrog(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
)
{
    /* Declare variables */
    int return_code;

    real *restrict x = system->x;
    real *restrict v = system->v;
    const int objects_count = system->objects_count;

    const real dt = integrator_param->dt;
    const int64 n_steps = simulation_param->n_steps_;

    const int storing_freq = storing_param->storing_freq;

    /* Allocate memory */
    real *restrict temp_x = malloc(objects_count * 3 * sizeof(real));
    real *restrict temp_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict a = malloc(objects_count * 3 * sizeof(real));
    
    // Compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    // Check if memory allocation is successful
    if (
        !temp_x 
        || !temp_v 
        || !a
        || !x_err_comp_sum 
        || !v_err_comp_sum
    )
    {
        return_code = ERROR_LEAPFROG_MEMORY_ALLOC;
        goto err_memory;
    }

    // Compute initial acceleration and v_1/2
    return_code = acceleration(
        a,
        system,
        acceleration_param
    );
    if (return_code != SUCCESS)
    {
        goto err_acceleration;
    }

    memcpy(temp_v, v, objects_count * 3 * sizeof(real));
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            v_err_comp_sum[i * 3 + j] += 0.5 * a[i * 3 + j] * dt;
            v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
            v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];
        }
    }

    /* Main Loop */
    for (int64 count = 1; count <= n_steps; count++)
    {
        // Calculate x_1
        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                x[i * 3 + j] = temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += temp_x[i * 3 + j] - x[i * 3 + j];
            }
        }

        // Calculate v_1+1/2
        return_code = acceleration(
            a,
            system,
            acceleration_param
        );
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;
                v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];
            }
        }

        /* Update time */
        *(simulation_status->t) = count * dt;

        /* Store solution */
        if (count % storing_freq == 0)
        {
            // Get v_1 from v_1+1/2
            memcpy(temp_v, v, objects_count * 3 * sizeof(real));
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    v[i * 3 + j] -= 0.5 * a[i * 3 + j] * dt;
                }
            }
            return_code = store_solution_step(
                storing_param,
                system,
                simulation_status,
                solutions
            );
            if (return_code != SUCCESS)
            {
                goto err_store_solution;
            }

            // Restore v_1+1/2
            memcpy(v, temp_v, objects_count * 3 * sizeof(real));
        }

        /* Check user interrupt */
        if (*(settings->is_exit))
        {
            return_code = ERROR_USER_INTERRUPT;
            goto err_user_interrupt;
        }
    }

    return SUCCESS;

err_user_interrupt:
err_store_solution:
err_acceleration:
err_memory:
    free(temp_x);
    free(temp_v);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    return return_code;
}
