#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

/**
 * \brief Compute the velocity kick
 * 
 * \param objects_count Number of objects in the system
 * \param jacobi_v Array of Jacobi velocity vectors
 * \param a Array of acceleration vectors
 * \param dt Time step of the system
 */
void whfast_kick(
    int objects_count,
    real *restrict jacobi_v,
    real *restrict a,
    real dt
);

/** 
 * \brief Compute the position drift
 * 
 * \param objects_count Number of objects in the system
 * \param jacobi_x Array of Jacobi position vectors
 * \param jacobi_v Array of Jacobi velocity vectors
 * \param m Array of masses
 * \param eta Array of cumulative masses
 * \param G Gravitational constant
 * \param dt Time step of the system
 * \param kepler_tol Tolerance for solving Kepler's equation
 * \param kepler_max_iter Maximum number of iterations in solving Kepler's equation
 * \param kepler_auto_remove Flag to indicate whether to remove objects
 *                           that failed to converge in Kepler's equation
 * \param kepler_failed_bool_array Array of flags to indicate whether 
 *                                 an object failed to converge in Kepler's equation
 * \param kepler_failed_flag Flag to indicate whether any object failed to converge
 *                           in Kepler's equation
 */
void whfast_drift(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    const real *restrict m,
    const real *restrict eta,
    real G,
    real dt,
    real kepler_tol,
    int kepler_max_iter,
    bool kepler_auto_remove,
    real kepler_auto_remove_tol,
    bool *restrict kepler_failed_bool_array,
    bool *restrict kepler_failed_flag
);

/**
 * \brief Compute the pairwise gravitational acceleration
 * 
 * \details This is a brute-force pairwise calculation
 *          of gravitational acceleration between all objects,
 *          which is O(n^2) complexity.
 * 
 * \param objects_count Number of objects in the system
 * \param jacobi_x Array of Jacobi position vectors
 * \param x Array of position vectors
 * \param a Array of acceleration vectors
 * \param m Array of masses
 * \param eta Array of cumulative masses
 * \param G Gravitational constant
 */
void whfast_acceleration_pairwise(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    const real *restrict eta,
    real G
);

/**
 * \brief Compute the gravitational acceleration, separating massive and massless objects
 * 
 * \details This function calculates the gravitational acceleration
 *          between massive and massless objects separately.
 *          This is an O(m^2 + mn) complexity calculation,
 *          where m and n are the number of massive and massless 
 *          objects, respectively.
 * 
 * \param objects_count Number of objects in the system
 * \param jacobi_x Array of Jacobi position vectors
 * \param x Array of position vectors
 * \param a Array of acceleration vectors
 * \param m Array of masses
 * \param eta Array of cumulative masses
 * \param G Gravitational constant
 */
void whfast_acceleration_massless(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    const real *restrict eta,
    real G
);

/**
 * \brief Transform Cartesian coordinates to Jacobi coordinates
 * 
 * \param objects_count Number of objects in the system
 * \param jacobi_x Array of Jacobi position vectors to be stored
 * \param jacobi_v Array of Jacobi velocity vectors to be stored
 * \param x Array of position vectors
 * \param v Array of velocity vectors
 * \param m Array of masses
 * \param eta Array of cumulative masses
 */
void cartesian_to_jacobi(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real *restrict x,
    real *restrict v,
    const real *restrict m,
    const real *restrict eta
);

/**
 * \brief Transform Jacobi coordinates to Cartesian coordinates
 * 
 * \param objects_count Number of objects in the system
 * \param jacobi_x Array of Jacobi position vectors
 * \param jacobi_v Array of Jacobi velocity vectors
 * \param x Array of position vectors to be stored
 * \param v Array of velocity vectors to be stored
 * \param m Array of masses
 * \param eta Array of cumulative masses
 */
void jacobi_to_cartesian(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real *restrict x,
    real *restrict v,
    const real *restrict m,
    const real *restrict eta
);

/**
 * \brief Compute the Stumpff functions c0, c1, c2, and c3 for a given argument z
 * 
 * \param z Input value
 * \param c0 Pointer to store c0
 * \param c1 Pointer to store c1
 * \param c2 Pointer to store c2
 * \param c3 Pointer to store c3
 */
void stumpff_functions(
    real z,
    real *restrict c0,
    real *restrict c1,
    real *restrict c2,
    real *restrict c3
);

/**
 * \brief Propagate the position and velocity vectors using Kepler's equation
 * 
 * \param i Index of the object
 * \param jacobi_x Array of Jacobi position vectors
 * \param jacobi_v Array of Jacobi velocity vectors
 * \param gm Product of gravitational constant and total mass between
 *           the object and the central object
 * \param dt Time step of the system
 * \param kepler_tol Tolerance for solving Kepler's equation
 * \param kepler_max_iter Maximum number of iterations in solving Kepler's equation
 * \param kepler_auto_remove Flag to indicate whether to remove objects
 *                           that failed to converge in Kepler's equation
 * \param kepler_failed_bool_array Array of flags to indicate whether
 *                                an object failed to converge in Kepler's equation
 * \param kepler_failed_flag Flag to indicate whether any object failed to converge
 *                           in Kepler's equation
 */
void propagate_kepler(
    int i,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real gm,
    real dt,
    real kepler_tol,
    int kepler_max_iter,
    bool kepler_auto_remove,
    real kepler_auto_remove_tol,
    bool *restrict kepler_failed_bool_array,
    bool *restrict kepler_failed_flag
);

/**
 * \brief WHFast integrator
 * 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors
 * \param v Array of velocity vectors
 * \param m Array of masses
 * \param G Gravitational constant
 * \param dt Time step of the system
 * \param acceleration_method Method to calculate acceleration
 * \param npts Number of time steps to be integrated
 * \param store_npts Number of points to be stored
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param kepler_tol Tolerance for solving Kepler's equation
 * \param kepler_max_iter Maximum number of iterations in solving Kepler's equation
 * \param kepler_auto_remove Flag to indicate whether to remove objects
 *                           that failed to converge in Kepler's equation
 * \param kepler_actual_objects_count Pointer to the actual number of objects
 *                                    after clearing objects
 * \param flush Flag to indicate whether to store solution into data file directly
 * \param flush_path Path to the file to store the solution
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to flag that indicates whether user sent 
 *                KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
WIN32DLL_API int whfast(
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
    real kepler_tol,
    int kepler_max_iter,
    bool kepler_auto_remove,
    real kepler_auto_remove_tol,
    int *restrict kepler_actual_objects_count,
    const bool flush,
    const char *restrict flush_path,
    Solutions *restrict solution,
    bool *restrict is_exit
)
{   
    void (*whfast_acceleration)(
        int objects_count,
        real *restrict jacobi_x,
        real *restrict x,
        real *restrict a,
        const real *restrict m,
        const real *restrict eta,
        real G
    );

    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        whfast_acceleration = whfast_acceleration_pairwise;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        whfast_acceleration = whfast_acceleration_massless;
    }
    else
    {
        printf("Error: acceleration method not recognized\n");
        goto err_acc_method;
    }

    // Allocate memory for calculation
    real *restrict jacobi_x = calloc(objects_count * 3, sizeof(real));
    real *restrict jacobi_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict temp_jacobi_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict a = malloc(objects_count * 3 * sizeof(real));
    real *restrict eta = malloc(objects_count * sizeof(real));

    if (
        !jacobi_x ||
        !jacobi_v ||
        !temp_jacobi_v ||
        !a ||
        !eta
    )
    {
        printf("Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Auto remove objects that failed to converge in Kepler's equation
    bool *kepler_failed_bool_array = NULL;
    bool kepler_failed_flag = false;
    if (kepler_auto_remove)
    {
        kepler_failed_bool_array = calloc(objects_count, sizeof(bool));

        if (!kepler_failed_bool_array)
        {
            printf("Error: Failed to allocate memory for kepler_failed_bool_array\n");
            goto err_kepler_memory;
        }
    }

    // Allocate memory for solution output
    FILE *flush_file = NULL;
    double *sol_state = NULL;
    double *sol_time = NULL;
    double *sol_dt = NULL;
    if (flush)
    {
        flush_file = fopen(flush_path, "w");

        if (!flush_file)
        {
            printf("Error: Failed to open file for flushing\n");
            goto err_flush_file;
        }

        // Initial value
        write_to_csv_file(flush_file, 0.0, dt, objects_count, x, v, m, G);
    } 
    else
    {
        sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
        sol_time = malloc(store_npts * sizeof(double));
        sol_dt = malloc(store_npts * sizeof(double));

        if (!sol_state || !sol_time || !sol_dt)
        {
            printf("Error: Failed to allocate memory for solution output\n");
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

    eta[0] = m[0];
    for (int i = 1; i < objects_count; i++)
    {
        eta[i] = eta[i - 1] + m[i];
    }
    cartesian_to_jacobi(objects_count, jacobi_x, jacobi_v, x, v, m, eta);
    whfast_acceleration(objects_count, jacobi_x, x, a, m, eta, G);
    whfast_kick(objects_count, jacobi_v, a, 0.5 * dt);

    // Main Loop
    for (int64 count = 1; count <= npts; count++)
    {
        whfast_drift(
            objects_count,
            jacobi_x,
            jacobi_v,
            m,
            eta,
            G,
            dt,
            kepler_tol,
            kepler_max_iter,
            kepler_auto_remove,
            kepler_auto_remove_tol,
            kepler_failed_bool_array,
            &kepler_failed_flag
        );
        // Remove objects that failed to converge in Kepler's equation
        //
        // IMPORTANT:
        // It is important to remove objects right after drift step
        // for large N, otherwise one hyperbolic object could pollute
        // the data of the whole system with nan values
        if (kepler_auto_remove && kepler_failed_flag)
        {
            kepler_failed_flag = false;

            // Remove objects
            int kepler_remove_count = 0;
            for (int i = 1; i < objects_count; i++)
            {
                if (kepler_failed_bool_array[i])
                {
                    kepler_failed_bool_array[i] = false;
                    kepler_remove_count++;
                    printf("kepler_auto_remove: Object %d with mass %f removed\n", i, m[i]);
                }
                else if (kepler_remove_count > 0)
                {
                    memcpy(&jacobi_x[(i - kepler_remove_count) * 3], &jacobi_x[i * 3], 3 * sizeof(real));
                    memcpy(&jacobi_v[(i - kepler_remove_count) * 3], &jacobi_v[i * 3], 3 * sizeof(real));
                    m[i - kepler_remove_count] = m[i];
                }
            }
    
            objects_count -= kepler_remove_count;
            for (int i = 1; i < objects_count; i++)
            {
                eta[i] = eta[i - 1] + m[i];
            }
            printf("%d object(s) removed in total. Remaining objects: %d\n", kepler_remove_count, objects_count);
        }
        jacobi_to_cartesian(objects_count, jacobi_x, jacobi_v, x, v, m, eta);
        whfast_acceleration(objects_count, jacobi_x, x, a, m, eta, G);
        whfast_kick(objects_count, jacobi_v, a, dt);

        // Store solution
        if (count % store_every_n == 0)
        {
            memcpy(temp_jacobi_v, jacobi_v, objects_count * 3 * sizeof(real));
            whfast_kick(objects_count, temp_jacobi_v, a, -0.5 * dt);
            jacobi_to_cartesian(objects_count, jacobi_x, temp_jacobi_v, x, v, m, eta);
            
            if (flush)
            {
                write_to_csv_file(flush_file, dt * count, dt, objects_count, x, v, m, G);
            }
            else
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
    *kepler_actual_objects_count = objects_count;

    free(jacobi_x);
    free(jacobi_v);
    free(temp_jacobi_v);
    free(a);
    free(eta);
    if (kepler_auto_remove)
    {
        free(kepler_failed_bool_array);
    }

    if (flush)
    {
        fclose(flush_file);
    }
    else
    {
        solution->sol_state = sol_state;
        solution->sol_time = sol_time;
        solution->sol_dt = sol_dt;
    }

    return 0;

err_user_exit: // User sends KeyboardInterrupt in main thread
err_flush_file:
err_sol_output_memory:
    if (flush)
    {
        fclose(flush_file);
    }
    else
    {
        free(sol_state);
        free(sol_time);
        free(sol_dt);
    }
err_kepler_memory:
    if (kepler_auto_remove)
    {
        free(kepler_failed_bool_array);
    }
err_calc_memory:
    free(jacobi_x);
    free(jacobi_v);
    free(temp_jacobi_v);
    free(a);
    free(eta);
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

void whfast_kick(
    int objects_count,
    real *restrict jacobi_v,
    real *restrict a,
    real dt
)
{
    for (int i = 0; i < objects_count; i++)
    {
        jacobi_v[i * 3 + 0] += a[i * 3 + 0] * dt;
        jacobi_v[i * 3 + 1] += a[i * 3 + 1] * dt;
        jacobi_v[i * 3 + 2] += a[i * 3 + 2] * dt;
    }
}

void whfast_drift(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    const real *restrict m,
    const real *restrict eta,
    real G,
    real dt,
    real kepler_tol,
    int kepler_max_iter,
    bool kepler_auto_remove,
    real kepler_auto_remove_tol,
    bool *restrict kepler_failed_bool_array,
    bool *restrict kepler_failed_flag
)
{
    for (int i = 1; i < objects_count; i++)
    {
        real gm = G * m[0] * eta[i] / eta[i - 1];
        propagate_kepler(
            i,
            jacobi_x,
            jacobi_v,
            gm,
            dt,
            kepler_tol,
            kepler_max_iter,
            kepler_auto_remove,
            kepler_auto_remove_tol,
            kepler_failed_bool_array,
            kepler_failed_flag
        );
    }
}

void whfast_acceleration_pairwise(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    const real *restrict eta,
    real G
)
{
    real aux[3];
    real temp_vec[3];
    real temp_vec_norm;
    real temp_vec_norm_cube;
    real temp_jacobi_norm;
    real temp_jacobi_norm_cube;
    for (int i = 1; i < objects_count; i++)
    {
        // Calculate x_0i
        temp_vec[0] = x[i * 3 + 0] - x[0];
        temp_vec[1] = x[i * 3 + 1] - x[1];
        temp_vec[2] = x[i * 3 + 2] - x[2];

        temp_vec_norm = vec_norm(temp_vec, 3);
        temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;
        temp_jacobi_norm = vec_norm(&jacobi_x[i * 3], 3);
        temp_jacobi_norm_cube = temp_jacobi_norm * temp_jacobi_norm * temp_jacobi_norm;
        for (int j = 0; j < 3; j++)
        {
            a[i * 3 + j] = G * m[0] * eta[i] / eta[i - 1]
            * (
                jacobi_x[i * 3 + j] / temp_jacobi_norm_cube
                - temp_vec[j] / temp_vec_norm_cube
            );
        }

        for (int j = 1; j < i; j++)
        {
            // Calculate x_ji
            temp_vec[0] = x[i * 3 + 0] - x[j * 3 + 0];
            temp_vec[1] = x[i * 3 + 1] - x[j * 3 + 1];
            temp_vec[2] = x[i * 3 + 2] - x[j * 3 + 2];

            temp_vec_norm = vec_norm(temp_vec, 3);
            temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

            aux[0] += G * m[j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[i * 3 + 0] -= aux[0] * eta[i] / eta[i - 1];
        a[i * 3 + 1] -= aux[1] * eta[i] / eta[i - 1];
        a[i * 3 + 2] -= aux[2] * eta[i] / eta[i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = i + 1; j < objects_count; j++)
        {
            // Calculate x_ij
            temp_vec[0] = x[j * 3 + 0] - x[i * 3 + 0];
            temp_vec[1] = x[j * 3 + 1] - x[i * 3 + 1];
            temp_vec[2] = x[j * 3 + 2] - x[i * 3 + 2];

            temp_vec_norm = vec_norm(temp_vec, 3);
            temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

            aux[0] += G * m[j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[i * 3 + 0] += aux[0];
        a[i * 3 + 1] += aux[1];
        a[i * 3 + 2] += aux[2];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = 0; j < i; j++)
        {
            for (int k = i + 1; k < objects_count; k++)
            {
                // Calculate x_jk
                temp_vec[0] = x[k * 3 + 0] - x[j * 3 + 0];
                temp_vec[1] = x[k * 3 + 1] - x[j * 3 + 1];
                temp_vec[2] = x[k * 3 + 2] - x[j * 3 + 2];

                temp_vec_norm = vec_norm(temp_vec, 3);
                temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

                aux[0] += G * m[j] * m[k] * temp_vec[0] / temp_vec_norm_cube;
                aux[1] += G * m[j] * m[k] * temp_vec[1] / temp_vec_norm_cube;
                aux[2] += G * m[j] * m[k] * temp_vec[2] / temp_vec_norm_cube;
            }
        }
        a[i * 3 + 0] -= aux[0] / eta[i - 1];
        a[i * 3 + 1] -= aux[1] / eta[i - 1];
        a[i * 3 + 2] -= aux[2] / eta[i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;
    }
}

void whfast_acceleration_massless(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    const real *restrict eta,
    real G
)
{
    real aux[3];
    real temp_vec[3];
    real temp_vec_norm;
    real temp_vec_norm_cube;
    real temp_jacobi_norm;
    real temp_jacobi_norm_cube;

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

    int idx_i;
    int idx_j;
    int idx_k;

    // Acceleration calculation for massive objects
    for (int i = 1; i < massive_objects_count; i++)
    {
        idx_i = massive_indices[i];

        // Calculate x_0i
        temp_vec[0] = x[idx_i * 3 + 0] - x[0];
        temp_vec[1] = x[idx_i * 3 + 1] - x[1];
        temp_vec[2] = x[idx_i * 3 + 2] - x[2];

        temp_vec_norm = vec_norm(temp_vec, 3);
        temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;
        temp_jacobi_norm = vec_norm(&jacobi_x[idx_i * 3], 3);
        temp_jacobi_norm_cube = temp_jacobi_norm * temp_jacobi_norm * temp_jacobi_norm;
        for (int j = 0; j < 3; j++)
        {
            a[idx_i * 3 + j] = G * m[0] * eta[idx_i] / eta[idx_i - 1]
            * (
                jacobi_x[idx_i * 3 + j] / temp_jacobi_norm_cube
                - temp_vec[j] / temp_vec_norm_cube
            );
        }

        for (int j = 1; j < i; j++)
        {
            idx_j = massive_indices[j];

            // Calculate x_ji
            temp_vec[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            temp_vec[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            temp_vec[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];

            temp_vec_norm = vec_norm(temp_vec, 3);
            temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

            aux[0] += G * m[idx_j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[idx_j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[idx_j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[idx_i * 3 + 0] -= aux[0] * eta[idx_i] / eta[idx_i - 1];
        a[idx_i * 3 + 1] -= aux[1] * eta[idx_i] / eta[idx_i - 1];
        a[idx_i * 3 + 2] -= aux[2] * eta[idx_i] / eta[idx_i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = i + 1; j < massive_objects_count; j++)
        {
            idx_j = massive_indices[j];

            // Calculate x_ij
            temp_vec[0] = x[idx_j * 3 + 0] - x[idx_i * 3 + 0];
            temp_vec[1] = x[idx_j * 3 + 1] - x[idx_i * 3 + 1];
            temp_vec[2] = x[idx_j * 3 + 2] - x[idx_i * 3 + 2];

            temp_vec_norm = vec_norm(temp_vec, 3);
            temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

            aux[0] += G * m[idx_j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[idx_j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[idx_j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[idx_i * 3 + 0] += aux[0];
        a[idx_i * 3 + 1] += aux[1];
        a[idx_i * 3 + 2] += aux[2];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = 0; j < i; j++)
        {
            idx_j = massive_indices[j];

            for (int k = i + 1; k < massive_objects_count; k++)
            {
                idx_k = massive_indices[k];

                // Calculate x_jk
                temp_vec[0] = x[idx_k * 3 + 0] - x[idx_j * 3 + 0];
                temp_vec[1] = x[idx_k * 3 + 1] - x[idx_j * 3 + 1];
                temp_vec[2] = x[idx_k * 3 + 2] - x[idx_j * 3 + 2];

                temp_vec_norm = vec_norm(temp_vec, 3);
                temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

                aux[0] += G * m[idx_j] * m[idx_k] * temp_vec[0] / temp_vec_norm_cube;
                aux[1] += G * m[idx_j] * m[idx_k] * temp_vec[1] / temp_vec_norm_cube;
                aux[2] += G * m[idx_j] * m[idx_k] * temp_vec[2] / temp_vec_norm_cube;
            }
        }
        a[idx_i * 3 + 0] -= aux[0] / eta[idx_i - 1];
        a[idx_i * 3 + 1] -= aux[1] / eta[idx_i - 1];
        a[idx_i * 3 + 2] -= aux[2] / eta[idx_i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;
    }

    // Acceleration calculation for massless objects
    for (int i = 0; i < massless_objects_count; i++)
    {
        idx_i = massless_indices[i];
        if (idx_i == 0)
        {
            continue;
        }
        
        // Calculate x_0i
        temp_vec[0] = x[idx_i * 3 + 0] - x[0];
        temp_vec[1] = x[idx_i * 3 + 1] - x[1];
        temp_vec[2] = x[idx_i * 3 + 2] - x[2];

        temp_vec_norm = vec_norm(temp_vec, 3);
        temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;
        temp_jacobi_norm = vec_norm(&jacobi_x[idx_i * 3], 3);
        temp_jacobi_norm_cube = temp_jacobi_norm * temp_jacobi_norm * temp_jacobi_norm;
        for (int j = 0; j < 3; j++)
        {
            a[idx_i * 3 + j] = G * m[0]
            * (
                jacobi_x[idx_i * 3 + j] / temp_jacobi_norm_cube
                - temp_vec[j] / temp_vec_norm_cube
            );
        }

        for (int j = 1; j < massive_objects_count; j++)
        {
            idx_j = massive_indices[j];
            if (idx_j >= idx_i)
            {
                break;
            }

            // Calculate x_ji
            temp_vec[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            temp_vec[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            temp_vec[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];

            temp_vec_norm = vec_norm(temp_vec, 3);
            temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

            aux[0] += G * m[idx_j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[idx_j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[idx_j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[idx_i * 3 + 0] -= aux[0];
        a[idx_i * 3 + 1] -= aux[1];
        a[idx_i * 3 + 2] -= aux[2];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = 1; j < massive_objects_count; j++)
        {
            idx_j = massive_indices[j];
            if (idx_j <= idx_i)
            {
                continue;
            }

            // Calculate x_ij
            temp_vec[0] = x[idx_j * 3 + 0] - x[idx_i * 3 + 0];
            temp_vec[1] = x[idx_j * 3 + 1] - x[idx_i * 3 + 1];
            temp_vec[2] = x[idx_j * 3 + 2] - x[idx_i * 3 + 2];

            temp_vec_norm = vec_norm(temp_vec, 3);
            temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

            aux[0] += G * m[idx_j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[idx_j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[idx_j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[idx_i * 3 + 0] += aux[0];
        a[idx_i * 3 + 1] += aux[1];
        a[idx_i * 3 + 2] += aux[2];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = 0; j < massive_objects_count; j++)
        {
            idx_j = massive_indices[j];
            if (idx_j >= idx_i)
            {
                break;
            }

            for (int k = j + 1; k < massive_objects_count; k++)
            {
                idx_k = massive_indices[k];
                if (idx_k <= idx_i)
                {
                    continue;
                }

                // Calculate x_jk
                temp_vec[0] = x[idx_k * 3 + 0] - x[idx_j * 3 + 0];
                temp_vec[1] = x[idx_k * 3 + 1] - x[idx_j * 3 + 1];
                temp_vec[2] = x[idx_k * 3 + 2] - x[idx_j * 3 + 2];

                temp_vec_norm = vec_norm(temp_vec, 3);
                temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

                aux[0] += G * m[idx_j] * m[idx_k] * temp_vec[0] / temp_vec_norm_cube;
                aux[1] += G * m[idx_j] * m[idx_k] * temp_vec[1] / temp_vec_norm_cube;
                aux[2] += G * m[idx_j] * m[idx_k] * temp_vec[2] / temp_vec_norm_cube;
            }
        }
        a[idx_i * 3 + 0] -= aux[0] / eta[idx_i - 1];
        a[idx_i * 3 + 1] -= aux[1] / eta[idx_i - 1];
        a[idx_i * 3 + 2] -= aux[2] / eta[idx_i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;
    }

    free(massive_indices);
    free(massless_indices);
}

void cartesian_to_jacobi(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real *restrict x,
    real *restrict v,
    const real *restrict m,
    const real *restrict eta
)
{
    real x_cm[3];
    real v_cm[3];

    x_cm[0] = m[0] * x[0];
    x_cm[1] = m[0] * x[1];
    x_cm[2] = m[0] * x[2];

    v_cm[0] = m[0] * v[0];
    v_cm[1] = m[0] * v[1];
    v_cm[2] = m[0] * v[2];

    for (int i = 1; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            jacobi_x[i * 3 + j] = x[i * 3 + j] - x_cm[j] / eta[i - 1];
            jacobi_v[i * 3 + j] = v[i * 3 + j] - v_cm[j] / eta[i - 1];
        
            x_cm[j] = x_cm[j] * (1.0 + m[i] / eta[i - 1]) + m[i] * jacobi_x[i * 3 + j];
            v_cm[j] = v_cm[j] * (1.0 + m[i] / eta[i - 1]) + m[i] * jacobi_v[i * 3 + j];
        }
    }

    jacobi_x[0] = x_cm[0] / eta[objects_count - 1];
    jacobi_x[1] = x_cm[1] / eta[objects_count - 1];
    jacobi_x[2] = x_cm[2] / eta[objects_count - 1];

    jacobi_v[0] = v_cm[0] / eta[objects_count - 1];
    jacobi_v[1] = v_cm[1] / eta[objects_count - 1];
    jacobi_v[2] = v_cm[2] / eta[objects_count - 1];
}

void jacobi_to_cartesian(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real *restrict x,
    real *restrict v,
    const real *restrict m,
    const real *restrict eta
)
{
    real x_cm[3];
    real v_cm[3];

    x_cm[0] = eta[objects_count - 1] * jacobi_x[0];
    x_cm[1] = eta[objects_count - 1] * jacobi_x[1];
    x_cm[2] = eta[objects_count - 1] * jacobi_x[2];

    v_cm[0] = eta[objects_count - 1] * jacobi_v[0];
    v_cm[1] = eta[objects_count - 1] * jacobi_v[1];
    v_cm[2] = eta[objects_count - 1] * jacobi_v[2];

    for (int i = (objects_count - 1); i > 0; i--)
    {
        for (int j = 0; j < 3; j++)
        {
            x_cm[j] = (x_cm[j] - m[i] * jacobi_x[i * 3 + j]) / eta[i];
            v_cm[j] = (v_cm[j] - m[i] * jacobi_v[i * 3 + j]) / eta[i];

            x[i * 3 + j] = jacobi_x[i * 3 + j] + x_cm[j];
            v[i * 3 + j] = jacobi_v[i * 3 + j] + v_cm[j];

            x_cm[j] = eta[i - 1] * x_cm[j];
            v_cm[j] = eta[i - 1] * v_cm[j];
        }
    }
    
    x[0] = x_cm[0] / m[0];
    x[1] = x_cm[1] / m[0];
    x[2] = x_cm[2] / m[0];

    v[0] = v_cm[0] / m[0];
    v[1] = v_cm[1] / m[0];
    v[2] = v_cm[2] / m[0];
}

void stumpff_functions(
    real z,
    real *restrict c0,
    real *restrict c1,
    real *restrict c2,
    real *restrict c3
)
{
    // Reduce the argument
    int n = 0;
    while (fabs(z) > 0.1)
    {
        z /= 4.0;
        n++;
    }

    // Compute stumpff functions
    *c3 = (1.0 - z / 20.0 * (1.0 - z / 42.0 * (1.0 - z / 72.0 * (1.0 - z / 110.0 \
                 * (1.0 - z / 156.0 * (1.0 - z / 210.0)))))
             ) / 6.0;
    *c2 = (1.0 - z / 12.0 * (1.0 - z / 30.0 * (1.0 - z / 56.0 * (1.0 - z / 90.0 \
                * (1.0 - z / 132.0 * (1.0 - z / 182.0)))))
            ) / 2.0;
    *c1 = 1.0 - z * (*c3);
    *c0 = 1.0 - z * (*c2);

    // Half-angle formulae to recover the actual argument
    while (n > 0)
    {
        *c3 = ((*c2) + (*c0) * (*c3)) / 4.0;
        *c2 = ((*c1) * (*c1)) / 2.0;
        *c1 = (*c0) * (*c1);
        *c0 = 2.0 * (*c0) * (*c0) - 1.0;
        n--;
    }
}

void propagate_kepler(
    int i,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real gm,
    real dt,
    real kepler_tol,
    int kepler_max_iter,
    bool kepler_auto_remove,
    real kepler_auto_remove_tol,
    bool *restrict kepler_failed_bool_array,
    bool *restrict kepler_failed_flag
)
{
    real x[3];
    real v[3];
    memcpy(x, &jacobi_x[i * 3], 3 * sizeof(real));
    memcpy(v, &jacobi_v[i * 3], 3 * sizeof(real));

    real x_norm = vec_norm(x, 3);
    real v_norm = vec_norm(v, 3);

    // Radial velocity
    real radial_v = vec_dot(x, v, 3) / x_norm; 

    real alpha = 2.0 * gm / x_norm - (v_norm * v_norm);

    /* Solve Kepler's equation with Newton-Raphson method */
    
    // Initial guess
    real s = dt / x_norm;

    // Solve Kepler's equation
    real c0 = 0.0;
    real c1 = 0.0;
    real c2 = 0.0;
    real c3 = 0.0;
    int j = 0;
    for (; j < kepler_max_iter; j++)
    {
        // Compute Stumpff functions
        stumpff_functions(alpha * (s * s), &c0, &c1, &c2, &c3);

        // Evaluate Kepler's equation and its derivative
        real F = (
            x_norm * s * c1
            + x_norm * radial_v * (s * s) * c2
            + gm * (s * s * s) * c3
            - dt
        );
        real dF = (
            x_norm * c0
            + x_norm * radial_v * s * c1
            + gm * (s * s) * c2
        );

        // Advance step
        real ds = -F / dF;
        s += ds;

        // Check convergence
        if (fabs(ds) < kepler_tol)
        {
            break;
        }
    }

    // The raidal distance is equal to the derivative of F
    // real r = dF
    real r = x_norm * c0 + x_norm * radial_v * s * c1 + gm * (s * s) * c2;

    if (j == kepler_max_iter)
    {
        real error = (x_norm * s * c1
            + x_norm * radial_v * (s * s) * c2
            + gm * (s * s * s) * c3
            - dt) / r;
        printf("Warning: Kepler's equation did not converge\n");
        printf("Object index: %d, error = %23.15g\n", i, error);
        if (kepler_auto_remove && ((fabs(error) > kepler_auto_remove_tol) || isnan(error)))
        {
            kepler_failed_bool_array[i] = true;
            *kepler_failed_flag = true;
        }

        // Debug information
        // printf("Input: x: %23.15g %23.15g %23.15g, v: %23.15g %23.15g %23.15g, gm: %23.15g, dt: %23.15g\n", x[0], x[1], x[2], v[0], v[1], v[2], gm, dt);
    }

    // Evaluate f and g functions, together with their derivatives
    real f = 1.0 - gm * (s * s) * c2 / x_norm;
    real g = dt - gm * (s * s * s) * c3;

    real df = -gm * s * c1 / (r * x_norm);
    real dg = 1.0 - gm * (s * s) * c2 / r; 

    // Compute position and velocity vectors
    for (int j = 0; j < 3; j++)
    {
        jacobi_x[i * 3 + j] = f * x[j] + g * v[j];
        jacobi_v[i * 3 + j] = df * x[j] + dg * v[j];
    }
}