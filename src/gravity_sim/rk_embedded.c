#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"


WIN32DLL_API void rk_embedded(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    const real *restrict m, 
    real G, 
    real expected_time_scale,
    real *restrict t, 
    real *restrict dt,
    int power,
    int power_test,
    int len_coeff,
    const real (*restrict coeff)[len_coeff],
    int len_weights,
    const real *restrict weights,
    const real *restrict weights_test,
    int max_iteration,
    int min_iteration,
    real abs_tolerance,
    real rel_tolerance
)
{
    // Initialization
    real t0 = *t;
    int stages = len_weights;
    int min_power = fmin(power, power_test);

    real *error_estimation_delta_weights = malloc(len_weights * sizeof(real));
    for (int stage = 0; stage < stages; stage++)
    {
        error_estimation_delta_weights[stage] = weights[stage] - weights_test[stage];
    }

    // Safety factors for step-size control:
    real safety_fac_max = 6.0;
    real safety_fac_min = 0.33;
    real safety_fac = pow(0.38, (1.0 / (1.0 + (real) min_power)));

    // Initialize arrays and values
    real sum, error, dt_new; 
    real (*v_1)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*x_1)[3] = malloc(objects_count * 3 * sizeof(real));
    real *vk = malloc(stages * objects_count * 3 * sizeof(real));
    real *xk = malloc(stages * objects_count * 3 * sizeof(real));    
    real (*temp_a)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*error_estimation_delta_v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*error_estimation_delta_x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*tolerance_scale_v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*tolerance_scale_x)[3] = malloc(objects_count * 3 * sizeof(real));

    // Main Loop
    for (int i = 0; i < max_iteration; i++)
    {
        // Calculate xk and vk
        acceleration(objects_count, x, temp_a, m, G);
        memcpy(vk, temp_a, objects_count * 3 * sizeof(real));
        memcpy(xk, v, objects_count * 3 * sizeof(real));

        for (int stage = 1; stage < stages; stage++)
        {
            // Empty temp_v and temp_x
            for (int i = 0; i < objects_count; i++)
            {
                temp_v[i][0] = 0.0;
                temp_v[i][1] = 0.0;
                temp_v[i][2] = 0.0;
                temp_x[i][0] = 0.0;
                temp_x[i][1] = 0.0;
                temp_x[i][2] = 0.0;
            }           

            for (int j = 0; j < stage; j++)
            {
                for (int k = 0; k < objects_count; k++)
                {
                    temp_v[k][0] += coeff[stage - 1][j] * vk[j * objects_count * 3 + k * 3];
                    temp_v[k][1] += coeff[stage - 1][j] * vk[j * objects_count * 3 + k * 3 + 1];
                    temp_v[k][2] += coeff[stage - 1][j] * vk[j * objects_count * 3 + k * 3 + 2];
                    temp_x[k][0] += coeff[stage - 1][j] * xk[j * objects_count * 3 + k * 3];
                    temp_x[k][1] += coeff[stage - 1][j] * xk[j * objects_count * 3 + k * 3 + 1];
                    temp_x[k][2] += coeff[stage - 1][j] * xk[j * objects_count * 3 + k * 3 + 2];
                }
            }

            for (int k = 0; k < objects_count; k++)
            {
                temp_x[k][0] = x[k][0] + *dt * temp_x[k][0];
                temp_x[k][1] = x[k][1] + *dt * temp_x[k][1];
                temp_x[k][2] = x[k][2] + *dt * temp_x[k][2];
                temp_v[k][0] = v[k][0] + *dt * temp_v[k][0];
                temp_v[k][1] = v[k][1] + *dt * temp_v[k][1];
                temp_v[k][2] = v[k][2] + *dt * temp_v[k][2];
            }
            acceleration(objects_count, temp_x, temp_a, m, G);
            memcpy(&vk[stage * objects_count * 3], temp_a, objects_count * 3 * sizeof(real));
            memcpy(&xk[stage * objects_count * 3], temp_v, objects_count * 3 * sizeof(real));
        }

        // Empty temp_v, temp_x, error_estimation_delta_v, error_estimation_delta_x
        for (int i = 0; i < objects_count; i++)
        {
            temp_v[i][0] = 0.0;
            temp_v[i][1] = 0.0;
            temp_v[i][2] = 0.0;
            temp_x[i][0] = 0.0;
            temp_x[i][1] = 0.0;
            temp_x[i][2] = 0.0;
            error_estimation_delta_v[i][0] = 0.0;
            error_estimation_delta_v[i][1] = 0.0;
            error_estimation_delta_v[i][2] = 0.0;
            error_estimation_delta_x[i][0] = 0.0;
            error_estimation_delta_x[i][1] = 0.0;
            error_estimation_delta_x[i][2] = 0.0;
        }

        // Calculate x_1, v_1 and also delta x, delta v for error estimation
        for(int stage = 0; stage < stages; stage++)
        {
            for (int k = 0; k < objects_count; k++)
            {
                temp_v[k][0] += weights[stage] * vk[stage * objects_count * 3 + k * 3];
                temp_v[k][1] += weights[stage] * vk[stage * objects_count * 3 + k * 3 + 1];
                temp_v[k][2] += weights[stage] * vk[stage * objects_count * 3 + k * 3 + 2];
                temp_x[k][0] += weights[stage] * xk[stage * objects_count * 3 + k * 3];
                temp_x[k][1] += weights[stage] * xk[stage * objects_count * 3 + k * 3 + 1];
                temp_x[k][2] += weights[stage] * xk[stage * objects_count * 3 + k * 3 + 2];

                error_estimation_delta_v[k][0] += *dt * error_estimation_delta_weights[stage] * vk[stage * objects_count * 3 + k * 3];
                error_estimation_delta_v[k][1] += *dt * error_estimation_delta_weights[stage] * vk[stage * objects_count * 3 + k * 3 + 1];
                error_estimation_delta_v[k][2] += *dt * error_estimation_delta_weights[stage] * vk[stage * objects_count * 3 + k * 3 + 2];
                error_estimation_delta_x[k][0] += *dt * error_estimation_delta_weights[stage] * xk[stage * objects_count * 3 + k * 3];
                error_estimation_delta_x[k][1] += *dt * error_estimation_delta_weights[stage] * xk[stage * objects_count * 3 + k * 3 + 1];
                error_estimation_delta_x[k][2] += *dt * error_estimation_delta_weights[stage] * xk[stage * objects_count * 3 + k * 3 + 2];
            }
        }

        for (int k = 0; k < objects_count; k++)
        {
            v_1[k][0] = v[k][0] + *dt * temp_v[k][0];
            v_1[k][1] = v[k][1] + *dt * temp_v[k][1];
            v_1[k][2] = v[k][2] + *dt * temp_v[k][2];
            x_1[k][0] = x[k][0] + *dt * temp_x[k][0];
            x_1[k][1] = x[k][1] + *dt * temp_x[k][1];
            x_1[k][2] = x[k][2] + *dt * temp_x[k][2];
        }

        // Error calculation
        for (int k = 0; k < objects_count; k++)
        {
            tolerance_scale_v[k][0] = abs_tolerance + fmax(fabs(v[k][0]), fabs(v_1[k][0])) * rel_tolerance;
            tolerance_scale_v[k][1] = abs_tolerance + fmax(fabs(v[k][1]), fabs(v_1[k][1])) * rel_tolerance;
            tolerance_scale_v[k][2] = abs_tolerance + fmax(fabs(v[k][2]), fabs(v_1[k][2])) * rel_tolerance;
            tolerance_scale_x[k][0] = abs_tolerance + fmax(fabs(x[k][0]), fabs(x_1[k][0])) * rel_tolerance;
            tolerance_scale_x[k][1] = abs_tolerance + fmax(fabs(x[k][1]), fabs(x_1[k][1])) * rel_tolerance;
            tolerance_scale_x[k][2] = abs_tolerance + fmax(fabs(x[k][2]), fabs(x_1[k][2])) * rel_tolerance;
        }

        // Sum up all the elements of x/tol and v/tol, 
        // square and divide by the total number of elements
        sum = 0.0;
        for (int k = 0; k < objects_count; k++)
        {
            sum += (error_estimation_delta_v[k][0] / tolerance_scale_v[k][0]) * (error_estimation_delta_v[k][0] / tolerance_scale_v[k][0]);
            sum += (error_estimation_delta_v[k][1] / tolerance_scale_v[k][1]) * (error_estimation_delta_v[k][1] / tolerance_scale_v[k][1]);
            sum += (error_estimation_delta_v[k][2] / tolerance_scale_v[k][2]) * (error_estimation_delta_v[k][2] / tolerance_scale_v[k][2]); 
            sum += (error_estimation_delta_x[k][0] / tolerance_scale_x[k][0]) * (error_estimation_delta_x[k][0] / tolerance_scale_x[k][0]); 
            sum += (error_estimation_delta_x[k][1] / tolerance_scale_x[k][1]) * (error_estimation_delta_x[k][1] / tolerance_scale_x[k][1]); 
            sum += (error_estimation_delta_x[k][2] / tolerance_scale_x[k][2]) * (error_estimation_delta_x[k][2] / tolerance_scale_x[k][2]); 
        }
        error = sqrt(sum / (objects_count * 3 * 2));

        if (error <= 1 || *dt == expected_time_scale * 1e-12)
        {
            // Advance step
            *t += *dt; 
            memcpy(x, x_1, objects_count * 3 * sizeof(real));
            memcpy(v, v_1, objects_count * 3 * sizeof(real));
        }

        // Calculate dt
        if (error != 0.0)   // Prevent division by zero
        {
            dt_new = *dt * safety_fac / pow(error, (1.0 / (1.0 + (real) min_power)));
        }
        else
        {
            dt_new = *dt;
        }
        
        if (dt_new > safety_fac_max * *dt) 
        {
            *dt *= safety_fac_max;
        }
        else if (dt_new < safety_fac_min * *dt)
        {
            *dt *= safety_fac_min;
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
        if (i >= min_iteration && *t >= (t0 + expected_time_scale * 1e-5))
        {
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
            break;
        }
    }
}
