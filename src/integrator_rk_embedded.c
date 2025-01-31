/**
 * \file integrator_rk_embedded.c
 * \author Ching Yin Ng
 * \brief Function definitions for Embedded Runge-Kutta integrators
 * 
 * Function definitions for Embedded Runge-Kutta integrators. This is
 * a C implementation based on the reference:
 *   J. Roa, et al. Moving Planets Around: An Introduction to
 *   N-Body Simulations Applied to Exoplanetary Systems*, MIT
 *   Press, 2020
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

/**
 * \brief Get the order of the Embedded RK integrator
 * 
 * \param integrator Name of the integrator
 * \param order Pointer to the order of the integrator
 * 
 * \retval SUCCESS If successful
 * \retval ERROR_UNKNOWN_RK_EMBEDDED_METHOD If integrator is not recognized
 */
IN_FILE int get_rk_embedded_order(
    const char *integrator,
    int *restrict order
)
{
    if (strcmp(integrator, "rkf45") == 0)
    {
        *order = 45;
        return SUCCESS;
    }
    else if (strcmp(integrator, "dopri") == 0)
    {
        *order = 54;
        return SUCCESS;
    }
    else if (strcmp(integrator, "dverk") == 0)
    {
        *order = 65;
        return SUCCESS;
    }
    else if (strcmp(integrator, "rkf78") == 0)
    {
        *order = 78;
        return SUCCESS;
    }
    else
    {
        return ERROR_UNKNOWN_RK_EMBEDDED_METHOD;
    }
}

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
 * \retval SUCCESS If successful
 * \retval ERROR_RK_EMBEDDED_BUTCHER_TABLEAUS_UNKNOWN_ORDER
 *         If order is not one of 45 / 54 / 78 / 65
 * \retval ERROR_RK_EMBEDDED_BUTCHER_TABLEAUS_MEMORY_ALLOC 
 *         If failed to allocate memory for coeff, weights and weights_test
 */
IN_FILE int rk_embedded_butcher_tableaus(
    const int order,
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

    int return_code;
    *coeff = NULL;
    *weights = NULL;
    *weights_test = NULL;

    switch (order)
    {
        // RUNGE-KUTTA-FEHLBERG 4(5)
        case 45:
            // Order
            *power = 4;
            *power_test = 5;
            // nodes = np.array([1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5])
            *coeff = malloc(25 * sizeof(real));
            *len_weights = 6;
            *weights = malloc(*len_weights * sizeof(real));
            *weights_test = malloc(6 * sizeof(real));
            if (!*coeff || !*weights || !*weights_test)
            {
                goto error_memory;
            }
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
            memcpy(
                *weights,
                (real [6]) {
                    25.0L / 216.0L, 0.0L, 1408.0L / 2565.0L, 2197.0L / 4104.0L, -0.2L, 0.0L
                },
                6 * sizeof(real)
            );
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
            *len_weights = 7;
            *weights = malloc(*len_weights * sizeof(real));
            *weights_test = malloc(7 * sizeof(real));

            if (!*coeff || !*weights || !*weights_test)
            {
                goto error_memory;
            }

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
            memcpy(
                *weights,
                (real [7]) {
                    35.0L / 384.0L, 0.0L, 500.0L / 1113.0L, 125.0L / 192.0L, -2187.0L / 6784.0L, 11.0L / 84.0L, 0.0L
                },
                7 * sizeof(real)
            );
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
            *len_weights = 13;
            *weights = malloc(*len_weights * sizeof(real));
            *weights_test = malloc(13 * sizeof(real));
            if (!*coeff || !*weights || !*weights_test)
            {
                goto error_memory;
            }

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
            memcpy(
                *weights,
                (real [13]) {
                    41.0L / 840.0L, 0.0L, 0.0L, 0.0L, 0.0L, 34.0L / 105.0L, 9.0L / 35.0L, 9.0L / 35.0L, 9.0L / 280.0L, 9.0L / 280.0L, 41.0L / 840.0L, 0.0L, 0.0L
                },
                13 * sizeof(real)
            );
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
            *len_weights = 8;
            *weights = malloc(*len_weights * sizeof(real));
            *weights_test = malloc(8 * sizeof(real));
            if (!*coeff || !*weights || !*weights_test)
            {
                goto error_memory;
            }
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
            memcpy(
                *weights,
                (real [8]) {
                    3.0L / 40.0L, 0.0L, 875.0L / 2244.0L, 23.0L / 72.0L, 264.0L / 1955.0L, 0.0L, 125.0L / 11592.0L, 43.0L / 616.0L
                },
                8 * sizeof(real)
            );
            memcpy(
                *weights_test,
                (real [8]) {
                    13.0L / 160.0L, 0.0L, 2375.0L / 5984.0L, 5.0L / 16.0L, 12.0L / 85.0L, 3.0L / 44.0L, 0.0L, 0.0L
                },
                8 * sizeof(real)
            );

            break;

        default:
            return_code = ERROR_RK_EMBEDDED_BUTCHER_TABLEAUS_UNKNOWN_ORDER;
            goto error_code;
    }

    return SUCCESS;

error_memory:
    free(*coeff);
    free(*weights);
    free(*weights_test);
error_code:
    return return_code;
}

/**
 * \brief Calculate initial time step for Embedded Runge-Kutta integrator
 * 
 * \param rel_tolerance Relative tolerance
 * \param abs_tolerance Absolute tolerance
 * \param initial_dt Pointer to the initial time step
 * \param power Power of the integrator
 * \param system Pointer to the system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \note Modified to return dt * 1e-2 since this function gives initial dt thats too large
 * 
 * \retval SUCCESS If successful
 * \retval ERROR_RK_EMBEDDED_INITIAL_DT_MEMORY_ALLOC
 *         If failed to allocate memory for calculation
 * \retval ERROR_RK_EMBEDDED_INITIAL_DT_NEGATIVE
 *        If initial_dt is negative
 */
IN_FILE int rk_embedded_initial_dt(
    real rel_tolerance,
    real abs_tolerance,
    real *restrict initial_dt,
    const int power,
    System *system,
    AccelerationParam *acceleration_param
)
{
    int return_code;

    const int objects_count = system->objects_count;
    real *restrict x = system->x;
    real *restrict v = system->v;

    /* Allocate memory and declare variables */
    real *restrict tolerance_scale_x = malloc(objects_count * 3 * sizeof(real));
    real *restrict tolerance_scale_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict x_1 = malloc(objects_count * 3 * sizeof(real));
    real *restrict v_1 = malloc(objects_count * 3 * sizeof(real));
    real *restrict a_1 = malloc(objects_count * 3 * sizeof(real));
    real *restrict a = malloc(objects_count * 3 * sizeof(real));
    if (
        !tolerance_scale_x ||
        !tolerance_scale_v ||
        !x_1 ||
        !v_1 ||
        !a_1 ||
        !a
    )
    {
        return_code = ERROR_RK_EMBEDDED_INITIAL_DT_MEMORY_ALLOC;
        goto error_memory;
    }

    real sum_0 = 0.0;
    real sum_1 = 0.0;
    real sum_2 = 0.0;
    real d_0;
    real d_1;
    real d_2;
    real dt_0;
    real dt_1;
    real dt;
    System system_1 = {
        .objects_count = objects_count,
        .x = x_1,
        .v = v_1,
        .m = system->m,
        .G = system->G,
    };

    /* Compute acceleration */
    return_code = acceleration(
        a,
        system,
        acceleration_param
    );
    if (return_code != SUCCESS)
    {
        goto error_acc;
    }

    /**
     *  tolerance_scale_x = abs_tol + rel_tol * abs(x)
     *  tolerance_scale_v = abs_tol + rel_tol * abs(v)
     */
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            tolerance_scale_x[i * 3 + j] = abs_tolerance + rel_tolerance * fabs(x[i * 3 + j]);
            tolerance_scale_v[i * 3 + j] = abs_tolerance + rel_tolerance * fabs(v[i * 3 + j]);
        }
    }

    /**
     *  sum_0 = sum(square(x / tolerance_scale_x)) + sum(square(v / tolerance_scale_v))
     *  sum_1 = sum(square(v / tolerance_scale_x)) + sum(square(a / tolerance_scale_x))
     */
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

    d_0 = sqrt(sum_0 / (objects_count * 6));
    d_1 = sqrt(sum_1 / (objects_count * 6));

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
            x_1[i * 3 + j] = x[i * 3 + j] + (dt_0 / 100.0) * v[i * 3 + j];
            v_1[i * 3 + j] = v[i * 3 + j] + (dt_0 / 100.0) * a[i * 3 + j];
        }
    }

    return_code = acceleration(
        a_1,
        &system_1,
        acceleration_param
    );
    if (return_code != SUCCESS)
    {
        goto error_acc;
    }

    /* Calculate d_2 to measure how much the derivatives have changed */

    /**
     * sum_2 = sum(square((v_1 - v) / tolerance_scale_x)) + sum(square((a_1 - a) / tolerance_scale_v))
     */
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sum_2 += ((v_1[i * 3 + j] - v[i * 3 + j]) / tolerance_scale_x[i * 3 + j]) * ((v_1[i * 3 + j] - v[i * 3 + j]) / tolerance_scale_x[i * 3 + j]);
            sum_2 += ((a_1[i * 3 + j] - a[i * 3 + j]) / tolerance_scale_v[i * 3 + j]) * ((a_1[i * 3 + j] - a[i * 3 + j]) / tolerance_scale_v[i * 3 + j]);
        }
    }
    d_2 = sqrt(sum_2 / (objects_count * 6)) / dt_0;

    if (fmax(d_1, d_2) <= 1e-15)
    {
        dt_1 = fmax(1e-6L, dt_0 * 1e-3L);
    }
    {
        dt_1 = pow((0.01L / fmax(d_1, d_2)), (1.0L / (1 + power)));
    }
    dt = fmin(100.0L * dt_0, dt_1);

    /* the factor 1e-2 is a modification made by me */
    *initial_dt = dt * 1e-2;
    if (*initial_dt <= 0.0)
    {
        return_code = ERROR_RK_EMBEDDED_INITIAL_DT_NON_POSITIVE;
        goto error_dt;
    }

    free(tolerance_scale_x);
    free(tolerance_scale_v);
    free(x_1);
    free(v_1);
    free(a_1);
    free(a);

    return SUCCESS;

error_dt:
error_acc:
error_memory:
    free(tolerance_scale_x);
    free(tolerance_scale_v);
    free(x_1);
    free(v_1);
    free(a_1);
    free(a);
    return return_code;
}

WIN32DLL_API int rk_embedded(
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
    int return_code;

    /* Initialization */
    int order;
    int power;
    int power_test;
    real *coeff = NULL;
    int len_weights;
    real *weights = NULL;
    real *weights_test = NULL;

    return_code = get_rk_embedded_order(
        integrator_param->integrator,
        &order
    );
    if (return_code != SUCCESS)
    {
        goto error_order;
    }

    return_code = rk_embedded_butcher_tableaus(
        order,
        &power,
        &power_test,
        &coeff,
        &len_weights,
        &weights,
        &weights_test
    );
    if (return_code != SUCCESS)
    {
        goto error_butcher_tableaus;
    }

    int stages = len_weights;
    int min_power = power < power_test ? power : power_test;

    real *restrict error_estimation_delta_weights = malloc(len_weights * sizeof(real));
    if (!error_estimation_delta_weights)
    {
        return_code = ERROR_RK_EMBEDDED_MEMORY_ALLOC;
        goto error_error_estimation_delta_weights_memory_alloc;
    }

    for (int stage = 0; stage < stages; stage++)
    {
        error_estimation_delta_weights[stage] = weights[stage] - weights_test[stage];
    }

    /* tolerance */
    real abs_tolerance = integrator_param->tolerance;
    real rel_tolerance = integrator_param->tolerance;

    /* Safety factors for step-size control */
    real safety_fac_max = 6.0;
    real safety_fac_min = 0.33;
    real safety_fac = pow(0.38, (1.0 / (1.0 + (real) min_power)));

    /* Allocate memory and declare variables */
    real *restrict x = system->x;
    real *restrict v = system->v;
    real *restrict m = system->m;
    const real G = system->G;
    const int objects_count = system->objects_count;

    real dt;
    real *t = simulation_status->t;
    real tf = simulation_param->tf;

    real sum;
    real error;
    real dt_new;

    System temp_system = {
        .objects_count = objects_count,
        .x = NULL,
        .v = NULL,
        .m = m,
        .G = G,
    };

    real *restrict v_1 = malloc(objects_count * 3 * sizeof(real));
    real *restrict x_1 = malloc(objects_count * 3 * sizeof(real));
    real *restrict vk = malloc(stages * objects_count * 3 * sizeof(real));
    real *restrict xk = malloc(stages * objects_count * 3 * sizeof(real));    
    real *restrict temp_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict temp_x = malloc(objects_count * 3 * sizeof(real));
    real *restrict error_estimation_delta_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict error_estimation_delta_x = malloc(objects_count * 3 * sizeof(real));
    real *restrict tolerance_scale_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict tolerance_scale_x = malloc(objects_count * 3 * sizeof(real));

    // Compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict temp_x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict temp_v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (
        !v_1 || 
        !x_1 || 
        !vk || 
        !xk ||  
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
        return_code = ERROR_RK_EMBEDDED_MEMORY_ALLOC;
        goto error_memory;
    }

    /* Get initial dt */
    if (integrator_param->initial_dt > 0.0)
    {
        dt = integrator_param->initial_dt;
    }
    else
    {
        return_code = rk_embedded_initial_dt(
            rel_tolerance,
            abs_tolerance,
            &dt,
            power,
            system,
            acceleration_param
        );
        if (return_code != SUCCESS)
        {
            goto error_initial_dt;
        }
    }
    simulation_status->dt = dt;

    int64 count = 1; // Count for storing solutions, 1 for t_0
    int storing_freq = storing_param->storing_freq;

    while (*t < tf)
    {
        /* Compute xk and vk */
        return_code = acceleration(
            vk,
            system,
            acceleration_param
        );
        if (return_code != SUCCESS)
        {
            goto acc_error;
        }
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

            temp_system.x = temp_x;
            temp_system.v = temp_v;
            return_code = acceleration(
                &vk[stage * objects_count * 3],
                &temp_system,
                acceleration_param
            );
            if (return_code != SUCCESS)
            {
                goto acc_error;
            }
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

        /* Error calculation */
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

        /* Advance step */
        if (error <= 1 || dt <= tf * 1e-12)
        {
            *t += dt; 
            memcpy(x, x_1, objects_count * 3 * sizeof(real));
            memcpy(v, v_1, objects_count * 3 * sizeof(real));

            memcpy(x_err_comp_sum, temp_x_err_comp_sum, objects_count * 3 * sizeof(real));
            memcpy(v_err_comp_sum, temp_v_err_comp_sum, objects_count * 3 * sizeof(real));

            /* Store solution */
            if (count % storing_freq == 0)
            {
                if ((*solutions->sol_size_ + 1) > storing_param->max_sol_size_)
                {
                    return_code = extend_sol_memory_buffer(
                        solutions,
                        storing_param,
                        objects_count
                    );
                    if (return_code != SUCCESS)
                    {
                        goto err_store_solution;
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
            }
        }

        /* Calculate dt for next step */
        if (error < 1e-10)
        {
            error = 1e-10;  // Prevent error from being too small
        }
        dt_new = dt * safety_fac / pow(error, (1.0 / (1.0 + (real) min_power)));
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
        if (((*t) < tf) && (((*t) + dt) > tf))
        {
            dt = tf - (*t);
        }
        simulation_status->dt = dt;

        /* Check user interrupt */
        if (*(settings->is_exit))
        {
            return_code = ERROR_USER_INTERRUPT;
            goto err_user_interrupt;
        }
    }

    /* free memory */
    free(v_1);
    free(x_1);
    free(vk);
    free(xk);
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
    free(error_estimation_delta_weights);
    free(coeff);
    free(weights);
    free(weights_test); 

    return SUCCESS;

err_user_interrupt:
err_store_solution:
acc_error:
error_initial_dt:
error_memory:
    free(v_1);
    free(x_1);
    free(vk);
    free(xk);
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
error_error_estimation_delta_weights_memory_alloc:
    free(error_estimation_delta_weights);
error_butcher_tableaus:
error_order:
    return return_code;
}