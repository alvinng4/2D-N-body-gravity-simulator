#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

/**
 * \brief Embedded RK integrator
 * 
 * \param order Order of the integrator
 * \param system Name of the system
 * \param t Pointer to the current simulation time
 * \param tf Total time to be integrated
 * \param input_abs_tolerance Absolute tolerance of the integrator
 * \param input_rel_tolerance Relative tolerance of the integrator
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param custom_sys_x Array of the position vectors of customized system 
 * \param custom_sys_v Array of the velocity vectors of customized system 
 * \param custom_sys_m Array of the masses of customized system 
 * \param custom_sys_G Gravitational constant of customized system
 * \param custom_sys_object_count Number of objects in customized system
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to determine whether user sent KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
int rk_embedded(
    int order,
    const char *restrict system,
    double *restrict t,
    double tf, 
    double input_abs_tolerance,
    double input_rel_tolerance,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count,
    Solutions *restrict solution,
    int *restrict is_exit
);

/**
 * \brief Butcher tableaus for Embedded RK integrator
 * 
 * \param order Order of the integrator, must be one of 45 / 54 / 78 / 65
 * \param power Power of the integrator
 * \param power_test Power for error calculation
 * \param coeff Pointer to array of coefficients for the integrator
 * \param len_weights Length of the weights array
 * \param weights Pointer to the array of weights for RK integrator
 * \param weights_test Pointer to the array of weights for error calculation
 * 
 * \retval 0 If successful
 * \retval 1 If failed to allocate memory for coeff, weights and weights_test
 */
int rk_embedded_butcher_tableaus(
    int order,
    int *restrict power,
    int *restrict power_test,
    real **coeff,
    int *restrict len_weights,
    real **weights,
    real **weights_test
);

/**
 * \brief Calculate the initial time step for Embedded RK integrator
 * 
 * \param objects_count Number of objects in the system
 * \param power Power of the integrator
 * \param x Array of position vectors of all objects
 * \param v Array of velocity vectors of all objects
 * \param a Array of acceleration vectors of all objects
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param abs_tolerance Absolute tolerance of the integrator
 * \param rel_tolerance Relative tolerance of the integrator
 * 
 * \warning: Modified to return dt * 1e-2 since this function gives initial dt thats too large
 * 
 * \return initial dt for Embedded RK integrator
 * \retval -1.0 if failed to allocate memory for calculation
 */
real rk_embedded_initial_dt(
    int objects_count,
    int power,
    real *restrict x,
    real *restrict v,
    real *restrict a,
    real *restrict m,
    real G,
    real abs_tolerance,
    real rel_tolerance
);


WIN32DLL_API int rk_embedded(
    int order,
    const char *restrict system,
    double *restrict t,
    double tf, 
    double input_abs_tolerance,
    double input_rel_tolerance,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count,
    Solutions *restrict solution,
    int *restrict is_exit
)
{   
    // Initialize default system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int intitial_sys_flag = initialize_system(system, &x, &v, &m, &objects_count, &G);

    // Custom system
    if (intitial_sys_flag == 2)
    {
        printf("Error: Failed to allocate memory to initialize default system\n");
        goto err_default_sys_memory;
    }
    else if (intitial_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

        if (!x || !v || !m)
        {
            printf("Error: Failed to allocate memory for customized system\n");
            goto err_custom_sys_memory;
        }

        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = custom_sys_x[i * 3 + j];
                v[i * 3 + j] = custom_sys_v[i * 3 + j];
            }
            m[i] = custom_sys_m[i];
        }
        G = custom_sys_G;
    }
    
    // Initialization
    int power;
    int power_test;
    real *coeff = NULL;
    int len_weights;
    real *weights = NULL;
    real *weights_test = NULL;
    int butcher_tableaus_flag = rk_embedded_butcher_tableaus(order, &power, &power_test, &coeff, &len_weights, &weights, &weights_test);
    if (butcher_tableaus_flag == 1)
    {
        goto err_rk_embedded_butcher_tableaus;
    }

    int stages = len_weights;
    int min_power = power < power_test ? power : power_test;

    real *error_estimation_delta_weights = malloc(len_weights * sizeof(real));
    if (!error_estimation_delta_weights)
    {
        printf("Error: Failed to allocate memory for error_estimation_delta_weights variable in rk_embedded()\n");
        goto err_error_estimation_delta_weights;
    }

    for (int stage = 0; stage < stages; stage++)
    {
        error_estimation_delta_weights[stage] = weights[stage] - weights_test[stage];
    }

    // Convert the tolerance from double to type real
    real abs_tolerance = input_abs_tolerance;
    real rel_tolerance = input_rel_tolerance;

    // Safety factors for step-size control:
    real safety_fac_max = 6.0;
    real safety_fac_min = 0.33;
    real safety_fac = pow(0.38, (1.0 / (1.0 + (real) min_power)));

    // Initialize memory for calculation
    real sum, error, dt_new; 
    real *a = malloc(objects_count * 3 * sizeof(real));
    real *v_1 = calloc(objects_count * 3, sizeof(real));
    real *x_1 = calloc(objects_count * 3, sizeof(real));
    real *vk = calloc(stages * objects_count * 3, sizeof(real));
    real *xk = calloc(stages * objects_count * 3, sizeof(real));    
    real *temp_a = calloc(objects_count * 3, sizeof(real));
    real *temp_v = calloc(objects_count * 3, sizeof(real));
    real *temp_x = calloc(objects_count * 3, sizeof(real));
    real *error_estimation_delta_v = calloc(objects_count * 3, sizeof(real));
    real *error_estimation_delta_x = calloc(objects_count * 3, sizeof(real));
    real *tolerance_scale_v = calloc(objects_count * 3, sizeof(real));
    real *tolerance_scale_x = calloc(objects_count * 3, sizeof(real));

    // Arrays for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *temp_x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *temp_v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (
        !a ||
        !v_1 || 
        !x_1 || 
        !vk || 
        !xk || 
        !temp_a || 
        !temp_v || 
        !temp_x || 
        !error_estimation_delta_v || 
        !error_estimation_delta_x || 
        !tolerance_scale_v || 
        !tolerance_scale_x || 
        !x_err_comp_sum || 
        !v_err_comp_sum || 
        !temp_x_err_comp_sum || 
        !temp_v_err_comp_sum
    ) 
    {
        printf("Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(NPTS * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(NPTS * sizeof(double));
    double *sol_dt = malloc(NPTS * sizeof(double));
    int buffer_size = NPTS;

    if (!sol_state || !sol_time || !sol_dt)
    {
        printf("Error: Failed to allocate memory for solution output\n");
        goto err_sol_output_memory;
    }
    
    // For realloc solution output
    double *temp_sol_state = NULL;
    double *temp_sol_time = NULL;
    double *temp_sol_dt = NULL;

    // Initial value
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sol_state[i * 3 + j] = x[i * 3 + j];
            sol_state[objects_count * 3 + i * 3 + j] = v[i * 3 + j];
        }
    }
    sol_time[0] = 0.0;
    real dt = rk_embedded_initial_dt(objects_count, power, x, v, a, m, G, abs_tolerance, rel_tolerance);
    if (dt == -1.0)
    {
        goto err_initial_dt_memory;
    }
    sol_dt[0] = dt;

    // Main Loop
    while (1)
    {
        // Calculate xk and vk
        acceleration(objects_count, x, vk, m, G);
        memcpy(xk, v, objects_count * 3 * sizeof(real));
        for (int stage = 1; stage < stages; stage++)
        {
            // Empty temp_v and temp_x
            for (int i = 0; i < objects_count; i++)
            {
                temp_v[i * 3 + 0] = 0.0;
                temp_v[i * 3 + 1] = 0.0;
                temp_v[i * 3 + 2] = 0.0;
                temp_x[i * 3 + 0] = 0.0;
                temp_x[i * 3 + 1] = 0.0;
                temp_x[i * 3 + 2] = 0.0;
            }       

            for (int i = 0; i < stage; i++)
            {
                for (int j = 0; j < objects_count; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        temp_v[j * 3 + k] += coeff[(stage - 1) * (stages - 1) + i] * vk[i * objects_count * 3 + j * 3 + k];
                        temp_x[j * 3 + k] += coeff[(stage - 1) * (stages - 1) + i] * xk[i * objects_count * 3 + j * 3 + k];
                    }
                }
            }

            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    temp_v[i * 3 + j] = v[i * 3 + j] + dt * temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                    temp_x[i * 3 + j] = x[i * 3 + j] + dt * temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                }
            }
            acceleration(objects_count, temp_x, &vk[stage * objects_count * 3], m, G);
            memcpy(&xk[stage * objects_count * 3], temp_v, objects_count * 3 * sizeof(real));
        }

        // Empty temp_v, temp_x, error_estimation_delta_v, error_estimation_delta_x
        for (int i = 0; i < objects_count; i++)
        {
            temp_v[i * 3 + 0] = 0.0;
            temp_v[i * 3 + 1] = 0.0;
            temp_v[i * 3 + 2] = 0.0;
            temp_x[i * 3 + 0] = 0.0;
            temp_x[i * 3 + 1] = 0.0;
            temp_x[i * 3 + 2] = 0.0;
            error_estimation_delta_v[i * 3 + 0] = 0.0;
            error_estimation_delta_v[i * 3 + 1] = 0.0;
            error_estimation_delta_v[i * 3 + 2] = 0.0;
            error_estimation_delta_x[i * 3 + 0] = 0.0;
            error_estimation_delta_x[i * 3 + 1] = 0.0;
            error_estimation_delta_x[i * 3 + 2] = 0.0;
        }

        // Calculate x_1, v_1 and also delta x, delta v for error estimation
        for(int stage = 0; stage < stages; stage++)
        {
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    temp_v[i * 3 + j] += weights[stage] * vk[stage * objects_count * 3 + i * 3 + j];
                    temp_x[i * 3 + j] += weights[stage] * xk[stage * objects_count * 3 + i * 3 + j];

                    error_estimation_delta_v[i * 3 + j] += dt * error_estimation_delta_weights[stage] * vk[stage * objects_count * 3 + i * 3 + j];
                    error_estimation_delta_x[i * 3 + j] += dt * error_estimation_delta_weights[stage] * xk[stage * objects_count * 3 + i * 3 + j];
                }
            }
        }

        memcpy(temp_x_err_comp_sum, x_err_comp_sum, objects_count * 3 * sizeof(real));
        memcpy(temp_v_err_comp_sum, v_err_comp_sum, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                temp_v_err_comp_sum[i * 3 + j] += dt * temp_v[i * 3 + j];
                temp_x_err_comp_sum[i * 3 + j] += dt * temp_x[i * 3 + j];

                v_1[i * 3 + j] = v[i * 3 + j] + temp_v_err_comp_sum[i * 3 + j];
                x_1[i * 3 + j] = x[i * 3 + j] + temp_x_err_comp_sum[i * 3 + j];

                temp_v_err_comp_sum[i * 3 + j] += v[i * 3 + j] - v_1[i * 3 + j];
                temp_x_err_comp_sum[i * 3 + j] += x[i * 3 + j] - x_1[i * 3 + j];
            }
        }

        // Error calculation
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                tolerance_scale_v[i * 3 + j] = abs_tolerance + fmax(fabs(v[i * 3 + j]), fabs(v_1[i * 3 + j])) * rel_tolerance;
                tolerance_scale_x[i * 3 + j] = abs_tolerance + fmax(fabs(x[i * 3 + j]), fabs(x_1[i * 3 + j])) * rel_tolerance;
            }
        }

        // Sum up all the elements of x/tol and v/tol, 
        // square and divide by the total number of elements
        sum = 0.0;
        for (int i = 0; i < objects_count; i++)
        {
            real temp;
            for (int j = 0; j < 3; j++)
            {
                temp = error_estimation_delta_v[i * 3 + j] / tolerance_scale_v[i * 3 + j];
                sum += temp * temp;
                temp = error_estimation_delta_x[i * 3 + j] / tolerance_scale_x[i * 3 + j];
                sum += temp * temp;
            }
        }
        error = sqrt(sum / (objects_count * 3 * 2));

        if (error <= 1 || dt == tf * 1e-12)
        {
            // Advance step
            *t += dt; 
            memcpy(x, x_1, objects_count * 3 * sizeof(real));
            memcpy(v, v_1, objects_count * 3 * sizeof(real));
            count += 1;

            memcpy(x_err_comp_sum, temp_x_err_comp_sum, objects_count * 3 * sizeof(real));
            memcpy(v_err_comp_sum, temp_v_err_comp_sum, objects_count * 3 * sizeof(real));

            // Store step
            if ((count + 1) % store_every_n == 0)
            {
                sol_time[*store_count + 1] = *t;
                sol_dt[*store_count + 1] = dt;
                memcpy(&sol_state[(*store_count + 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
                memcpy(&sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
                *store_count += 1;
            }

            else if (*t >= tf)
            {
                sol_time[*store_count] = *t;
                sol_dt[*store_count] = dt;
                memcpy(&sol_state[(*store_count + 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
                memcpy(&sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
            }

            // Check if user sends KeyboardInterrupt in main thread
            if (*is_exit)
            {
                goto err_user_exit;
            }

            // End simulation as t >= tf
            else if (*t >= tf)
            {
                free(x);
                free(v);
                free(a);
                free(coeff);
                free(weights);
                free(weights_test);
                free(error_estimation_delta_weights);
                free(v_1);
                free(x_1);
                free(vk);
                free(xk);
                free(temp_a);
                free(temp_v);
                free(temp_x);
                free(error_estimation_delta_v);
                free(error_estimation_delta_x);
                free(tolerance_scale_v);
                free(tolerance_scale_x);
                free(x_err_comp_sum);
                free(v_err_comp_sum);
                free(temp_x_err_comp_sum);
                free(temp_v_err_comp_sum);

                solution->sol_state = sol_state;
                solution->sol_time = sol_time;
                solution->sol_dt = sol_dt;
                solution->m = m;
                solution->G = G;
                solution->objects_count = objects_count;

                return 0;
            }

            // Check buffer size and extend if full
            if ((*store_count + 1) == buffer_size)
            {   
                buffer_size *= 2;
                temp_sol_state = realloc(sol_state, buffer_size * objects_count * 6 * sizeof(real));
                temp_sol_time = realloc(sol_time, buffer_size * sizeof(real));
                temp_sol_dt = realloc(sol_dt, buffer_size * sizeof(real));

                if (!temp_sol_state || !temp_sol_time || !temp_sol_dt)
                {
                    printf("Error: Failed to allocate extra memory to extend array for solution output\n");
                    goto err_sol_output_memory;
                }
                
                sol_state = temp_sol_state;
                sol_time = temp_sol_time;
                sol_dt = temp_sol_dt;
            }
        }

        // Calculate dt
        if (error != 0.0)   // Prevent division by zero
        {
            dt_new = dt * safety_fac / pow(error, (1.0 / (1.0 + (real) min_power)));
        }
        else
        {
            dt_new = dt;
        }
        
        if (dt_new > safety_fac_max * dt) 
        {
            dt *= safety_fac_max;
        }
        else if (dt_new < safety_fac_min * dt)
        {
            dt *= safety_fac_min;
        }
        else
        {
            dt = dt_new;
        }

        if (dt_new / tf < 1e-12)
        {
            dt = tf * 1e-12;
        }

        // Correct overshooting
        if (*t + dt > tf)
        {
            dt = tf - *t;
        }
    }

err_user_exit: // User sends KeyboardInterrupt in main thread
err_initial_dt_memory:
err_sol_output_memory:
    free(sol_state);
    free(sol_time);
    free(sol_dt);
err_calc_memory:
    free(a);
    free(v_1);
    free(x_1);
    free(vk);
    free(xk);
    free(temp_a);
    free(temp_v);
    free(temp_x);
    free(error_estimation_delta_v);
    free(error_estimation_delta_x);
    free(tolerance_scale_v);
    free(tolerance_scale_x);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    free(temp_x_err_comp_sum);
    free(temp_v_err_comp_sum);
err_error_estimation_delta_weights:
    free(error_estimation_delta_weights);
err_rk_embedded_butcher_tableaus:
    free(coeff);
    free(weights);
    free(weights_test);
err_custom_sys_memory:
    free(x);
    free(v);
    free(m);
err_default_sys_memory:    
    if (*is_exit)
    {
        return 2;   // User sends KeyboardInterrupt in main thread
    }
    else
    {
        return 1;
    }
}

WIN32DLL_API int rk_embedded_butcher_tableaus(
    int order,
    int *restrict power,
    int *restrict power_test,
    real **coeff,
    int *restrict len_weights,
    real **weights,
    real **weights_test
)
{
    /*  
    *   Select integrator
    *   45) Runge-Kutta-Fehleberg 4(5)
    *   54) Dormand-Prince 5(4)
    *   78) Runge-Kutta-Fehlberg 7(8)
    *   65) Verner's method 6(5), DVERK
    */

    switch (order)
    {
        // RUNGE-KUTTA-FEHLBERG 4(5)
        case 45:
            // Order
            *power = 4;
            *power_test = 5;
            // nodes = np.array([1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5])
            *coeff = malloc(25 * sizeof(real));
            memcpy(
                *coeff,
                (real [25]) {
                    1.0L / 4.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    3.0L / 32.0L, 9.0L / 32.0L, 0.0L, 0.0L, 0.0L,
                    1932.0L / 2197.0L, -7200.0L / 2197.0L, 7296.0L / 2197.0L, 0.0L, 0.0L,
                    439.0L / 216.0L, -8.0L, 3680.0L / 513.0L, -845.0L / 4104.0L, 0.0L,
                    -8.0L / 27.0L, 2.0L, -3544.0L / 2565.0L, 1859.0L / 4104.0L, -11.0L / 40.0L
                },
                25 * sizeof(real)
            );

            *len_weights = 6;
            *weights = malloc(*len_weights * sizeof(real));
            memcpy(
                *weights,
                (real [6]) {
                    25.0L / 216.0L, 0.0L, 1408.0L / 2565.0L, 2197.0L / 4104.0L, -0.2L, 0.0L
                },
                6 * sizeof(real)
            );

            *weights_test = malloc(6 * sizeof(real));
            memcpy(
                *weights_test,
                (real [6]) {
                    16.0L / 135.0L, 0.0L, 6656.0L / 12825.0L, 28561.0L / 56430.0L, -9.0L / 50.0L, 2.0L / 55.0L
                },
                6 * sizeof(real)
            );

            break;

        // DORMAND-PRINCE 5(4)
        case 54:
            // order
            *power = 5;
            *power_test = 4;
            // nodes = np.array([1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0])
            *coeff = malloc(36 * sizeof(real));
            memcpy(
                *coeff,
                (real [36]) {
                    1.0L / 5.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    3.0L / 40.0L, 9.0L / 40.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    44.0L / 45.0L, -56.0L / 15.0L, 32.0L / 9.0L, 0.0L, 0.0L, 0.0L,
                    19372.0L / 6561.0L, -25360.0L / 2187.0L, 64448.0L / 6561.0L, -212.0L / 729.0L, 0.0L, 0.0L,
                    9017.0L / 3168.0L, -355.0L / 33.0L, 46732.0L / 5247.0L, 49.0L / 176.0L, -5103.0L / 18656.0L, 0.0L,
                    35.0L / 384.0L, 0.0L, 500.0L / 1113.0L, 125.0L / 192.0L, -2187.0L / 6784.0L, 11.0L / 84.0L
                },
                36 * sizeof(real)
            );

            *len_weights = 7;
            *weights = malloc(*len_weights * sizeof(real));
            memcpy(
                *weights,
                (real [7]) {
                    35.0L / 384.0L, 0.0L, 500.0L / 1113.0L, 125.0L / 192.0L, -2187.0L / 6784.0L, 11.0L / 84.0L, 0.0L
                },
                7 * sizeof(real)
            );

            *weights_test = malloc(7 * sizeof(real));
            memcpy(
                *weights_test,
                (real [7]) {
                    5179.0L / 57600.0L, 0.0L, 7571.0L / 16695.0L, 393.0L / 640.0L, -92097.0L / 339200.0L, 187.0L / 2100.0L, 1.0L / 40.0L
                },
                7 * sizeof(real)
            );

            break;

        // RUNGE-KUTTA-FEHLBERG 7(8)
        case 78:
            // Order
            *power = 7;
            *power_test = 8;
            // nodes = np.array(
            //     [
            //         2.0 / 27.0,
            //         1.0 / 9.0,
            //         1.0 / 6.0,
            //         5.0 / 12.0,
            //         1.0 / 2.0,
            //         5.0 / 6.0,
            //         1.0 / 6.0,
            //         2.0 / 3.0,
            //         1.0 / 3.0,
            //         1.0,
            //         0.0,
            //         1.0,
            //     ]
            // )
            
            *coeff = malloc(144 * sizeof(real));
            memcpy(
                *coeff,
                (real [144]) {
                    2.0L / 27.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    1.0L / 36.0L, 1.0L / 12.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    1.0L / 24.0L, 0.0L, 1.0L / 8.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    5.0L / 12.0L, 0.0L, -25.0L / 16.0L, 25.0L / 16.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    1.0L / 20.0L, 0.0L, 0.0L, 1.0L / 4.0L, 1.0L / 5.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    -25.0L / 108.0L, 0.0L, 0.0L, 125.0L / 108.0L, -65.0L / 27.0L, 125.0L / 54.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    31.0L / 300.0L, 0.0L, 0.0L, 0.0L, 61.0L / 225.0L, -2.0L / 9.0L, 13.0L / 900.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    2.0L, 0.0L, 0.0L, -53.0L / 6.0L, 704.0L / 45.0L, -107.0L / 9.0L, 67.0L / 90.0L, 3.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    -91.0L / 108.0L, 0.0L, 0.0L, 23.0L / 108.0L, -976.0L / 135.0L, 311.0L / 54.0L, -19.0L / 60.0L, 17.0L / 6.0L, -1.0L / 12.0L, 0.0L, 0.0L, 0.0L,
                    2383.0L / 4100.0L, 0.0L, 0.0L, -341.0L / 164.0L, 4496.0L / 1025.0L, -301.0L / 82.0L, 2133.0L / 4100.0L, 45.0L / 82.0L, 45.0L / 164.0L, 18.0L / 41.0L, 0.0L, 0.0L,
                    3.0L / 205.0L, 0.0L, 0.0L, 0.0L, 0.0L, -6.0L / 41.0L, -3.0L / 205.0L, -3.0L / 41.0L, 3.0L / 41.0L, 6.0L / 41.0L, 0.0L, 0.0L,
                    -1777.0L / 4100.0L, 0.0L, 0.0L, -341.0L / 164.0L, 4496.0L / 1025.0L, -289.0L / 82.0L, 2193.0L / 4100.0L, 51.0L / 82.0L, 33.0L / 164.0L, 19.0L / 41.0L, 0.0L, 1.0L
                },
                144 * sizeof(real)
            );

            *len_weights = 13;
            *weights = malloc(*len_weights * sizeof(real));
            memcpy(
                *weights,
                (real [13]) {
                    41.0L / 840.0L, 0.0L, 0.0L, 0.0L, 0.0L, 34.0L / 105.0L, 9.0L / 35.0L, 9.0L / 35.0L, 9.0L / 280.0L, 9.0L / 280.0L, 41.0L / 840.0L, 0.0L, 0.0L
                },
                13 * sizeof(real)
            );

            *weights_test = malloc(13 * sizeof(real));
            memcpy(
                *weights_test,
                (real [13]) {
                    0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 34.0L / 105.0L, 9.0L / 35.0L, 9.0L / 35.0L, 9.0L / 280.0L, 9.0L / 280.0L, 0.0L, 41.0L / 840.0L, 41.0L / 840.0L
                },
                13 * sizeof(real)
            );

            break;

        // VERNER 6(5) DVERK
        case 65:
            // Order
            *power = 6;
            *power_test = 7;
            /* nodes = np.array(
            *     [1.0 / 6.0, 4.0 / 15.0, 2.0 / 3.0, 5.0 / 6.0, 1.0, 1.0 / 15.0, 1.0]
            * )
            */
            *coeff = malloc(49 * sizeof(real));
            memcpy(
                *coeff,
                (real [49]) {
                    1.0L / 6.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    4.0L / 75.0L, 16.0L / 75.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    5.0L / 6.0L, -8.0L / 3.0L, 5.0L / 2.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    -165.0L / 64.0L, 55.0L / 6.0L, -425.0L / 64.0L, 85.0L / 96.0L, 0.0L, 0.0L, 0.0L,
                    12.0L / 5.0L, -8.0L, 4015.0L / 612.0L, -11.0L / 36.0L, 88.0L / 255.0L, 0.0L, 0.0L,
                    -8263.0L / 15000.0L, 124.0L / 75.0L, -643.0L / 680.0L, -81.0L / 250.0L, 2484.0L / 10625.0L, 0.0L, 0.0L,
                    3501.0L / 1720.0L, -300.0L / 43.0L, 297275.0L / 52632.0L, -319.0L / 2322.0L, 24068.0L / 84065.0L, 0.0L, 3850.0L / 26703.0L
                },
                49 * sizeof(real)
            );

            *len_weights = 8;
            *weights = malloc(*len_weights * sizeof(real));
            memcpy(
                *weights,
                (real [8]) {
                    3.0L / 40.0L, 0.0L, 875.0L / 2244.0L, 23.0L / 72.0L, 264.0L / 1955.0L, 0.0L, 125.0L / 11592.0L, 43.0L / 616.0L
                },
                8 * sizeof(real)
            );

            *weights_test = malloc(8 * sizeof(real));
            memcpy(
                *weights_test,
                (real [8]) {
                    13.0L / 160.0L, 0.0L, 2375.0L / 5984.0L, 5.0L / 16.0L, 12.0L / 85.0L, 3.0L / 44.0L, 0.0L, 0.0L
                },
                8 * sizeof(real)
            );

            break;

        default:
            printf("Error: Invalid order for Embedded RK integrator\n");
            exit(EXIT_FAILURE);
            break;
    }

    if (!*coeff || !*weights || !weights_test)
    {
        printf("Error: Failed to allocate memory for rk_embedded_butcher_tableaus()\n");
        free(*coeff);
        free(*weights);
        free(*weights_test);

        return 1;
    }
    return 0;
}

WIN32DLL_API real rk_embedded_initial_dt(
    int objects_count,
    int power,
    real *restrict x,
    real *restrict v,
    real *restrict a,
    real *restrict m,
    real G,
    real abs_tolerance,
    real rel_tolerance
)
{
    real *tolerance_scale_x = malloc(objects_count * 3 * sizeof(real));
    real *tolerance_scale_v = malloc(objects_count * 3 * sizeof(real));
    real sum_0 = 0;
    real sum_1 = 0;
    real sum_2 = 0;
    real d_0;
    real d_1;
    real d_2;
    real dt_0;
    real dt_1;
    real dt;
    real *x_1 = malloc(objects_count * 3 * sizeof(real));
    real *v_1 = malloc(objects_count * 3 * sizeof(real));
    real *a_1 = malloc(objects_count * 3 * sizeof(real));

    if (
        !tolerance_scale_x ||
        !tolerance_scale_v ||
        !x_1 ||
        !v_1 ||
        !a_1
    )
    {
        printf("Error: Failed to allocate memory for rk_embedded_initial_dt()\n");
        free(tolerance_scale_x);
        free(tolerance_scale_v);
        free(x_1);
        free(v_1);
        free(a_1);
        return -1.0;
    }

    acceleration(objects_count, x, a, m, G);

    // tolerance_scale_x = abs_tol + rel_tol * abs(x)
    // tolerance_scale_v = abs_tol + rel_tol * abs(v)
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            tolerance_scale_x[i * 3 + j] = abs_tolerance + rel_tolerance * fabs(x[i * 3 + j]);
            tolerance_scale_v[i * 3 + j] = abs_tolerance + rel_tolerance * fabs(v[i * 3 + j]);
        }
    }

    // sum_0 = sum(square(x / tolerance_scale_x)) + sum(square(v / tolerance_scale_v))
    // sum_1 = sum(square(v / tolerance_scale_x)) + sum(square(a / tolerance_scale_x))
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sum_0 += (x[i * 3 + j] / tolerance_scale_x[i * 3 + j]) * (x[i * 3 + j] / tolerance_scale_x[i * 3 + j]);
            sum_0 += (v[i * 3 + j] / tolerance_scale_v[i * 3 + j]) * (v[i * 3 + j] / tolerance_scale_v[i * 3 + j]);
            sum_1 += (v[i * 3 + j] / tolerance_scale_x[i * 3 + j]) * (v[i * 3 + j] / tolerance_scale_x[i * 3 + j]);
            sum_1 += (a[i * 3 + j] / tolerance_scale_v[i * 3 + j]) * (a[i * 3 + j] / tolerance_scale_v[i * 3 + j]);
        }
    }

    d_0 = pow(sum_0 / (objects_count * 6), 0.5);
    d_1 = pow(sum_1 / (objects_count * 6), 0.5);

    if (d_0 < 1e-5 || d_1 < 1e-5)
    {
        dt_0 = 1e-4;
    }
    else
    {
        dt_0 = d_0 / d_1;
    }

    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            x_1[i * 3 + j] = x[i * 3 + j] + (dt_0 / 100.0L) * v[i * 3 + j];
            v_1[i * 3 + j] = v[i * 3 + j] + (dt_0 / 100.0L) * a[i * 3 + j];
        }
    }
    
    acceleration(objects_count, x_1, a_1, m, G);

    // Calculate d_2 to measure how much the derivatives have changed

    // sum_2 = sum(square((v_1 - v) / tolerance_scale_x)) + sum(square((a_1 - a) / tolerance_scale_v))
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sum_2 += ((v_1[i * 3 + j] - v[i * 3 + j]) / tolerance_scale_x[i * 3 + j]) * ((v_1[i * 3 + j] - v[i * 3 + j]) / tolerance_scale_x[i * 3 + j]);
            sum_2 += ((a_1[i * 3 + j] - a[i * 3 + j]) / tolerance_scale_v[i * 3 + j]) * ((a_1[i * 3 + j] - a[i * 3 + j]) / tolerance_scale_v[i * 3 + j]);
        }
    }
    d_2 = pow(sum_2 / (objects_count * 6), 0.5) / dt_0;

    if (fmax(d_1, d_2) <= 1e-15)
    {
        dt_1 = fmax(1e-6L, dt_0 * 1e-3L);
    }
    {
        dt_1 = pow((0.01L / fmax(d_1, d_2)), (1.0L / (1 + power)));
    }
    dt = fmin(100.0L * dt_0, dt_1);

    free(tolerance_scale_x);
    free(tolerance_scale_v);
    free(x_1);
    free(v_1);
    free(a_1);

    return dt * 1e-2;
}