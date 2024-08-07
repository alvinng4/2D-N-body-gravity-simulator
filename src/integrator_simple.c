#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"


/**
 * \brief Euler integrator
 * 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param v Array of velocity vectors of all objects
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param dt Time step of the system
 * \param npts Number of time steps to be integrated
 * \param acceleration_method Method to calculate acceleration
 * \param store_npts Number of points to be stored
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param storing_method Integer flag to indicate method of storing solution
 * \param flush_path Path to the file to store the solution
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to flag that indicates whether user sent 
 *                KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
WIN32DLL_API int euler(
    int objects_count,
    real *restrict x,
    real *restrict v,
    real *restrict m,
    real G,
    double dt,
    const char *restrict acceleration_method,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const int storing_method,
    const char *restrict flush_path,
    Solutions *restrict solution,
    bool *restrict is_exit
)
{   
    void (*acceleration)(
        int objects_count,
        real *restrict x,
        real *restrict a,
        const real *restrict m,
        real G
    );
    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        acceleration = acceleration_pairwise;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        acceleration = acceleration_massless;
    }
    else if (strcmp(acceleration_method, "barnes-hut") == 0)
    {
        acceleration = acceleration_barnes_hut;
    }
    else
    {
        fprintf(stderr, "Error: acceleration method not recognized\n");
        goto err_acc_method;
    }

    // Allocate memory for calculation
    real *restrict temp_x = malloc(objects_count * 3 * sizeof(real));
    real *restrict temp_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict a = malloc(objects_count * 3 * sizeof(real));
    
    // Allocate memory for compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (!temp_x || !temp_v || !a || !x_err_comp_sum || !v_err_comp_sum)
    {
        fprintf(stderr, "Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    FILE *flush_file = NULL;
    double *sol_state = NULL;
    double *sol_time = NULL;
    double *sol_dt = NULL;
    if (storing_method == 1)
    {
        flush_file = fopen(flush_path, "w");

        if (!flush_file)
        {
            fprintf(stderr, "Error: Failed to open file for flushing\n");
            goto err_flush_file;
        }

        // Initial value
        write_to_csv_file(flush_file, 0.0, dt, objects_count, x, v, m, G);
    } 
    else if (storing_method == 0)
    {
        sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
        sol_time = malloc(store_npts * sizeof(double));
        sol_dt = malloc(store_npts * sizeof(double));

        if (!sol_state || !sol_time || !sol_dt)
        {
            fprintf(stderr, "Error: Failed to allocate memory for solution output\n");
            goto err_sol_output_memory;
        }

        // Initial value
        memcpy(&sol_state[0], x, objects_count * 3 * sizeof(double));
        memcpy(&sol_state[objects_count * 3], v, objects_count * 3 * sizeof(double));
        sol_time[0] = 0.0;
        for (int i = 0; i < store_npts; i++)
        {
            sol_dt[i] = dt;
        }
    }

    // Main Loop
    for (int64 count = 1; count <= npts; count++)
    {
        acceleration(objects_count, x, a, m, G);

        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));

        // Calculation
        for (int i = 0; i < objects_count; i++) 
        {
            for (int j = 0; j < 3; j++) 
            {
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;

                x[i * 3 + j] = temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];

                x_err_comp_sum[i * 3 + j] += temp_x[i * 3 + j] - x[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];
            }
        }

        // Store solution
        if (count % store_every_n == 0)
        {
            if (storing_method == 1)
            {
                write_to_csv_file(flush_file, dt * count, dt, objects_count, x, v, m, G);
            }
            else if (storing_method == 0)
            {
                memcpy(&sol_state[*store_count * objects_count * 6], x, objects_count * 3 * sizeof(double));
                memcpy(&sol_state[*store_count * objects_count * 6 + objects_count * 3], v, objects_count * 3 * sizeof(double));
                sol_time[*store_count] = dt * count;
            }
            (*store_count)++;
        }

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            goto err_user_exit;
        }
    }

    // Exit after simulation is finished
    free(a);
    free(temp_x);
    free(temp_v);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    if (storing_method == 1)
    {
        fclose(flush_file);
    }
    else if (storing_method == 0)
    {
        solution->sol_state = sol_state;
        solution->sol_time = sol_time;
        solution->sol_dt = sol_dt;
    }

    return 0;

err_user_exit: // User sends KeyboardInterrupt in main thread
err_flush_file:
err_sol_output_memory:
    if (storing_method == 1)
    {
        fclose(flush_file);
    }
    else if (storing_method == 0)
    {
        free(sol_state);
        free(sol_time);
        free(sol_dt);
    }    
err_calc_memory:
    free(temp_x);
    free(temp_v);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
err_acc_method:
    if (*is_exit)
    {
        return 2;   // User sends KeyboardInterrupt in main thread
    }
    else
    {
        return 1;
    }
}

/**
 * \brief Euler-cromer integrator
 * 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param v Array of velocity vectors of all objects
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param dt Time step of the system
 * \param acceleration_method Method to calculate acceleration
 * \param npts Number of time steps to be integrated
 * \param store_npts Number of points to be stored
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param storing_method Integer flag to indicate method of storing solution
 * \param flush_path Path to the file to store the solution
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to flag that indicates whether user sent 
 *                KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
WIN32DLL_API int euler_cromer(
    int objects_count,
    real *restrict x,
    real *restrict v,
    real *restrict m,
    real G,
    double dt,
    const char *restrict acceleration_method,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const int storing_method,
    const char *restrict flush_path,
    Solutions *restrict solution,
    bool *restrict is_exit
)
{   
    void (*acceleration)(
        int objects_count,
        real *restrict x,
        real *restrict a,
        const real *restrict m,
        real G
    );
    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        acceleration = acceleration_pairwise;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        acceleration = acceleration_massless;
    }
    else if (strcmp(acceleration_method, "barnes-hut") == 0)
    {
        acceleration = acceleration_barnes_hut;
    }
    else
    {
        fprintf(stderr, "Error: acceleration method not recognized\n");
        goto err_acc_method;
    }

    // Allocate memory for calculation
    real *restrict temp_x = malloc(objects_count * 3 * sizeof(real));
    real *restrict temp_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict a = malloc(objects_count * 3 * sizeof(real));
    
    // Allocate memory for compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    
    if (!temp_x || !temp_v || !a || !x_err_comp_sum || !v_err_comp_sum)
    {
        fprintf(stderr, "Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    FILE *flush_file = NULL;
    double *sol_state = NULL;
    double *sol_time = NULL;
    double *sol_dt = NULL;
    if (storing_method == 1)
    {
        flush_file = fopen(flush_path, "w");

        if (!flush_file)
        {
            fprintf(stderr, "Error: Failed to open file for flushing\n");
            goto err_flush_file;
        }

        // Initial value
        write_to_csv_file(flush_file, 0.0, dt, objects_count, x, v, m, G);
    } 
    else if (storing_method == 0)
    {
        sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
        sol_time = malloc(store_npts * sizeof(double));
        sol_dt = malloc(store_npts * sizeof(double));

        if (!sol_state || !sol_time || !sol_dt)
        {
            fprintf(stderr, "Error: Failed to allocate memory for solution output\n");
            goto err_sol_output_memory;
        }

        // Initial value
        memcpy(&sol_state[0], x, objects_count * 3 * sizeof(double));
        memcpy(&sol_state[objects_count * 3], v, objects_count * 3 * sizeof(double));
        sol_time[0] = 0.0;
        for (int i = 0; i < store_npts; i++)
        {
            sol_dt[i] = dt;
        }
    }

    // Main Loop
    for (int64 count = 1; count <= npts; count++)
    {
        acceleration(objects_count, x, a, m, G);

        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                // Calculate v
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;
                v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];

                // Calculate x
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                x[i * 3 + j] = temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += temp_x[i * 3 + j] - x[i * 3 + j];
            }
        }

        // Store solution
        if (count % store_every_n == 0)
        {
            if (storing_method == 1)
            {
                write_to_csv_file(flush_file, dt * count, dt, objects_count, x, v, m, G);
            }
            else if (storing_method == 0)
            {
                memcpy(&sol_state[*store_count * objects_count * 6], x, objects_count * 3 * sizeof(double));
                memcpy(&sol_state[*store_count * objects_count * 6 + objects_count * 3], v, objects_count * 3 * sizeof(double));
                sol_time[*store_count] = dt * count;
            }
            (*store_count)++;
        }

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            goto err_user_exit;
        }
    }

    // Exit after simulation is finished
    free(a);
    free(temp_x);
    free(temp_v);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    if (storing_method == 1)
    {
        fclose(flush_file);
    }
    else if (storing_method == 0)
    {
        solution->sol_state = sol_state;
        solution->sol_time = sol_time;
        solution->sol_dt = sol_dt;
    }

    return 0;

err_user_exit: // User sends KeyboardInterrupt in main thread
err_flush_file:
err_sol_output_memory:
    if (storing_method == 1)
    {
        fclose(flush_file);
    }
    else if (storing_method == 0)
    {
        free(sol_state);
        free(sol_time);
        free(sol_dt);
    }    
err_calc_memory:
    free(temp_x);
    free(temp_v);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
err_acc_method:
    if (*is_exit)
    {
        return 2;   // User sends KeyboardInterrupt in main thread
    }
    else
    {
        return 1;
    }
}

/**
 * \brief RK4 integrator
 * 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param v Array of velocity vectors of all objects
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param dt Time step of the system
 * \param acceleration_method Method to calculate acceleration
 * \param npts Number of time steps to be integrated
 * \param store_npts Number of points to be stored
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param storing_method Integer flag to indicate method of storing solution
 * \param flush_path Path to the file to store the solution
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to flag that indicates whether user sent 
 *                KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
WIN32DLL_API int rk4(
    int objects_count,
    real *restrict x,
    real *restrict v,
    real *restrict m,
    real G,
    double dt,
    const char *restrict acceleration_method,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const int storing_method,
    const char *restrict flush_path,
    Solutions *restrict solution,
    int *restrict is_exit
)
{
    void (*acceleration)(
        int objects_count,
        real *restrict x,
        real *restrict a,
        const real *restrict m,
        real G
    );
    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        acceleration = acceleration_pairwise;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        acceleration = acceleration_massless;
    }
    else if (strcmp(acceleration_method, "barnes-hut") == 0)
    {
        acceleration = acceleration_barnes_hut;
    }
    else
    {
        fprintf(stderr, "Error: acceleration method not recognized\n");
        goto err_acc_method;
    }

    // Allocate memory for calculation
    real *restrict temp_x = malloc(objects_count * 3 * sizeof(real));
    real *restrict temp_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict a = malloc(objects_count * 3 * sizeof(real));
    real *restrict vk1 = malloc(objects_count * 3 * sizeof(real));
    real *restrict vk2 = malloc(objects_count * 3 * sizeof(real));
    real *restrict vk3 = malloc(objects_count * 3 * sizeof(real));
    real *restrict vk4 = malloc(objects_count * 3 * sizeof(real));
    real *restrict xk1 = malloc(objects_count * 3 * sizeof(real));
    real *restrict xk2 = malloc(objects_count * 3 * sizeof(real));
    real *restrict xk3 = malloc(objects_count * 3 * sizeof(real));
    real *restrict xk4 = malloc(objects_count * 3 * sizeof(real));

    // Allocate memory for compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (!temp_x || !temp_v || !a || !vk1 || !vk2 || !vk3 || !vk4 || !xk1 || !xk2 || !xk3 || !xk4 || !x_err_comp_sum || !v_err_comp_sum)
    {
        fprintf(stderr, "Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    FILE *flush_file = NULL;
    double *sol_state = NULL;
    double *sol_time = NULL;
    double *sol_dt = NULL;
    if (storing_method == 1)
    {
        flush_file = fopen(flush_path, "w");

        if (!flush_file)
        {
            fprintf(stderr, "Error: Failed to open file for flushing\n");
            goto err_flush_file;
        }

        // Initial value
        write_to_csv_file(flush_file, 0.0, dt, objects_count, x, v, m, G);
    } 
    else if (storing_method == 0)
    {
        sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
        sol_time = malloc(store_npts * sizeof(double));
        sol_dt = malloc(store_npts * sizeof(double));

        if (!sol_state || !sol_time || !sol_dt)
        {
            fprintf(stderr, "Error: Failed to allocate memory for solution output\n");
            goto err_sol_output_memory;
        }

        // Initial value
        memcpy(&sol_state[0], x, objects_count * 3 * sizeof(double));
        memcpy(&sol_state[objects_count * 3], v, objects_count * 3 * sizeof(double));
        sol_time[0] = 0.0;
        for (int i = 0; i < store_npts; i++)
        {
            sol_dt[i] = dt;
        }
    }

    // Main Loop
    for (int64 count = 1; count <= npts; count++)
    {   
        acceleration(objects_count, x, vk1, m, G);
        memcpy(xk1, v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                temp_x[j * 3 + k] = x[j * 3 + k] + 0.5 * xk1[j * 3 + k] * dt;
                temp_v[j * 3 + k] = v[j * 3 + k] + 0.5 * vk1[j * 3 + k] * dt;
            }
        }
        acceleration(objects_count, temp_x, vk2, m, G);
        memcpy(xk2, temp_v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                temp_x[j * 3 + k] = x[j * 3 + k] + 0.5 * xk2[j * 3 + k] * dt;
                temp_v[j * 3 + k] = v[j * 3 + k] + 0.5 * vk2[j * 3 + k] * dt;
            }
        }
        acceleration(objects_count, temp_x, vk3, m, G);
        memcpy(xk3, temp_v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                temp_x[j * 3 + k] = x[j * 3 + k] + xk3[j * 3 + k] * dt;
                temp_v[j * 3 + k] = v[j * 3 + k] + vk3[j * 3 + k] * dt;
            }
        }
        acceleration(objects_count, temp_x, vk4, m, G);
        memcpy(xk4, temp_v, objects_count * 3 * sizeof(real));

        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        memcpy(temp_x, x, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            for (int k = 0; k < 3; k++)
            {
                v_err_comp_sum[j * 3 + k] += (vk1[j * 3 + k] + 2 * vk2[j * 3 + k] + 2 * vk3[j * 3 + k] + vk4[j * 3 + k]) * dt / 6.0;
                x_err_comp_sum[j * 3 + k] += (xk1[j * 3 + k] + 2 * xk2[j * 3 + k] + 2 * xk3[j * 3 + k] + xk4[j * 3 + k]) * dt / 6.0;

                v[j * 3 + k] = temp_v[j * 3 + k] + v_err_comp_sum[j * 3 + k];
                x[j * 3 + k] = temp_x[j * 3 + k] + x_err_comp_sum[j * 3 + k];

                v_err_comp_sum[j * 3 + k] += temp_v[j * 3 + k] - v[j * 3 + k];
                x_err_comp_sum[j * 3 + k] += temp_x[j * 3 + k] - x[j * 3 + k];
            }
        }

        // Store solution
        if (count % store_every_n == 0)
        {
            if (storing_method == 1)
            {
                write_to_csv_file(flush_file, dt * count, dt, objects_count, x, v, m, G);
            }
            else if (storing_method == 0)
            {
                memcpy(&sol_state[*store_count * objects_count * 6], x, objects_count * 3 * sizeof(double));
                memcpy(&sol_state[*store_count * objects_count * 6 + objects_count * 3], v, objects_count * 3 * sizeof(double));
                sol_time[*store_count] = dt * count;
            }
            (*store_count)++;
        }

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            goto err_user_exit;
        }
    }

    // Exit after simulation is finished
    free(a);
    free(temp_x);
    free(temp_v);
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

    if (storing_method == 1)
    {
        fclose(flush_file);
    }
    else if (storing_method == 0)
    {
        solution->sol_state = sol_state;
        solution->sol_time = sol_time;
        solution->sol_dt = sol_dt;
    }

    return 0;

err_user_exit: // User sends KeyboardInterrupt in main thread
err_flush_file:
err_sol_output_memory:
    if (storing_method == 1)
    {
        fclose(flush_file);
    }
    else if (storing_method == 0)
    {
        free(sol_state);
        free(sol_time);
        free(sol_dt);
    }    
err_calc_memory:
    free(temp_x);
    free(temp_v);
    free(a);
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
err_acc_method:
    if (*is_exit)
    {
        return 2;   // User sends KeyboardInterrupt in main thread
    }
    else
    {
        return 1;
    }
}

/**
 * \brief LeapFrog integrator
 * 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param v Array of velocity vectors of all objects
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param dt Time step of the system
 * \param acceleration_method Method to calculate acceleration
 * \param npts Number of time steps to be integrated
 * \param store_npts Number of points to be stored
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param storing_method Integer flag to indicate method of storing solution
 * \param flush_path Path to the file to store the solution
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to flag that indicates whether user sent 
 *                KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
WIN32DLL_API int leapfrog(
    int objects_count,
    real *restrict x,
    real *restrict v,
    real *restrict m,
    real G,
    double dt,
    const char *restrict acceleration_method,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const int storing_method,
    const char *restrict flush_path,
    Solutions *restrict solution,
    int *restrict is_exit
)
{   
    void (*acceleration)(
        int objects_count,
        real *restrict x,
        real *restrict a,
        const real *restrict m,
        real G
    );
    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        acceleration = acceleration_pairwise;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        acceleration = acceleration_massless;
    }
    else if (strcmp(acceleration_method, "barnes-hut") == 0)
    {
        acceleration = acceleration_barnes_hut;
    }
    else
    {
        fprintf(stderr, "Error: acceleration method not recognized\n");
        goto err_acc_method;
    }

    // Allocate memory for calculation
    real *restrict temp_x = malloc(objects_count * 3 * sizeof(real));
    real *restrict temp_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict a_0 = malloc(objects_count * 3 * sizeof(real));
    real *restrict a_1 = malloc(objects_count * 3 * sizeof(real));

    // Allocate memory for compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (!temp_x || !temp_v || !a_0 || !a_1 || !x_err_comp_sum || !v_err_comp_sum)
    {
        fprintf(stderr, "Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    FILE *flush_file = NULL;
    double *sol_state = NULL;
    double *sol_time = NULL;
    double *sol_dt = NULL;
    if (storing_method == 1)
    {
        flush_file = fopen(flush_path, "w");

        if (!flush_file)
        {
            fprintf(stderr, "Error: Failed to open file for flushing\n");
            goto err_flush_file;
        }

        // Initial value
        write_to_csv_file(flush_file, 0.0, dt, objects_count, x, v, m, G);
    } 
    else if (storing_method == 0)
    {
        sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
        sol_time = malloc(store_npts * sizeof(double));
        sol_dt = malloc(store_npts * sizeof(double));

        if (!sol_state || !sol_time || !sol_dt)
        {
            fprintf(stderr, "Error: Failed to allocate memory for solution output\n");
            goto err_sol_output_memory;
        }

        // Initial value
        memcpy(&sol_state[0], x, objects_count * 3 * sizeof(double));
        memcpy(&sol_state[objects_count * 3], v, objects_count * 3 * sizeof(double));
        sol_time[0] = 0.0;
        for (int i = 0; i < store_npts; i++)
        {
            sol_dt[i] = dt;
        }
    }

    // Main Loop
    acceleration(objects_count, x, a_1, m, G);
    for (int64 count = 1; count <= npts; count++)
    {       
        // Use a_1 from last iteration as a_0
        memcpy(a_0, a_1, objects_count * 3 * sizeof(real));

        // Calculate x
        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt + 0.5 * a_0[i * 3 + j] * dt * dt;
                x[i * 3 + j] = temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += temp_x[i * 3 + j] - x[i * 3 + j];
            }
        }

        // Calculate v
        acceleration(objects_count, x, a_1, m, G);
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                v_err_comp_sum[i * 3 + j] += 0.5 * (a_0[i * 3 + j] + a_1[i * 3 + j]) * dt;
                v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];
            }
        }

        // Store solution
        if (count % store_every_n == 0)
        {
            if (storing_method == 1)
            {
                write_to_csv_file(flush_file, dt * count, dt, objects_count, x, v, m, G);
            }
            else if (storing_method == 0)
            {
                memcpy(&sol_state[*store_count * objects_count * 6], x, objects_count * 3 * sizeof(double));
                memcpy(&sol_state[*store_count * objects_count * 6 + objects_count * 3], v, objects_count * 3 * sizeof(double));
                sol_time[*store_count] = dt * count;
            }
            (*store_count)++;
        }

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            goto err_user_exit;
        }
    }

    // Exit after simulation is finished
    free(temp_x);
    free(temp_v);
    free(a_0);
    free(a_1);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    if (storing_method == 1)
    {
        fclose(flush_file);
    }
    else if (storing_method == 0)
    {
        solution->sol_state = sol_state;
        solution->sol_time = sol_time;
        solution->sol_dt = sol_dt;
    }

    return 0;

err_user_exit: // User sends KeyboardInterrupt in main thread
err_flush_file:
err_sol_output_memory:
    if (storing_method == 1)
    {
        fclose(flush_file);
    }
    else if (storing_method == 0)
    {
        free(sol_state);
        free(sol_time);
        free(sol_dt);
    }    
err_calc_memory:
    free(temp_x);
    free(temp_v);
    free(a_0);
    free(a_1);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
err_acc_method:
    if (*is_exit)
    {
        return 2;   // User sends KeyboardInterrupt in main thread
    }
    else
    {
        return 1;
    }
}
