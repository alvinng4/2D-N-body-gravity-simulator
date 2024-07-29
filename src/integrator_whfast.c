#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

#define TOL_KEPLER 1e-12

void whfast_kick(
    int objects_count,
    real *restrict jacobi_v,
    real *restrict a,
    real dt
);

void whfast_drift(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    const real *restrict m,
    const real *restrict eta,
    real G,
    real dt
);

void whfast_acceleration_pairwise(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict x,
    real *restrict a,
    const real *restrict m,
    const real *restrict eta,
    real G
);

void cartesian_to_jacobi(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real *restrict x,
    real *restrict v,
    const real *restrict m,
    const real *restrict eta
);

void jacobi_to_cartesian(
    int objects_count,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real *restrict x,
    real *restrict v,
    const real *restrict m,
    const real *restrict eta
);

void stumpff_functions(
    real psi,
    real *restrict c2,
    real *restrict c3
);

void propagate_kepler(
    int i,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real gm,
    real dt
);

/**
 * \brief WHFast integrator
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
    // else if (strcmp(acceleration_method, "massless") == 0)
    // {
    //     whfast_acceleration = whfast_acceleration_massless;
    // }
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
    
    // Allocate memory for compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (
        !jacobi_x ||
        !jacobi_v ||
        !temp_jacobi_v ||
        !a ||
        !eta ||
        !x_err_comp_sum ||
        !v_err_comp_sum
    )
    {
        printf("Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
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
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                sol_state[i * 3 + j] = x[i * 3 + j];
                sol_state[objects_count * 3 + i * 3 + j] = v[i * 3 + j];
            }
        }
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
        whfast_drift(objects_count, jacobi_x, jacobi_v, m, eta, G, dt);
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
                memcpy(&sol_state[*store_count * objects_count * 6], x, objects_count * 6 * sizeof(double));
                memcpy(&sol_state[*store_count * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
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
    free(jacobi_x);
    free(jacobi_v);
    free(temp_jacobi_v);
    free(a);
    free(eta);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

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
err_calc_memory:
    free(jacobi_x);
    free(jacobi_v);
    free(temp_jacobi_v);
    free(a);
    free(eta);
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
    real dt
)
{
    for (int i = 1; i < objects_count; i++)
    {
        real gm = G * m[0] * eta[i] / eta[i - 1];
        propagate_kepler(i, jacobi_x, jacobi_v, gm, dt);
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

            for (int k = 0; k < 3; k++)
            {
                aux[k] += G * m[j] * temp_vec[k] / temp_vec_norm_cube;
            } 
        }
        for (int j = 0; j < 3; j++)
        {
            a[i * 3 + j] -= aux[j] * eta[i] / eta[i - 1];

            // Empty aux array
            aux[j] = 0.0;
        }

        for (int j = i + 1; j < objects_count; j++)
        {
            // Calculate x_ij
            temp_vec[0] = x[j * 3 + 0] - x[i * 3 + 0];
            temp_vec[1] = x[j * 3 + 1] - x[i * 3 + 1];
            temp_vec[2] = x[j * 3 + 2] - x[i * 3 + 2];

            temp_vec_norm = vec_norm(temp_vec, 3);
            temp_vec_norm_cube = temp_vec_norm * temp_vec_norm * temp_vec_norm;

            for (int k = 0; k < 3; k++)
            {
                aux[k] += G * m[j] * temp_vec[k] / temp_vec_norm_cube;
            }
        }
        for (int j = 0; j < 3; j++)
        {
            a[i * 3 + j] += aux[j];

            // Empty aux array
            aux[j] = 0.0;
        }

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

                for (int l = 0; l < 3; l++)
                {
                    aux[l] += G * m[j] * m[k] * temp_vec[l] / temp_vec_norm_cube;
                }
            }
        }
        for (int j = 0; j < 3; j++)
        {
            a[i * 3 + j] -= aux[j] / eta[i - 1];

            // Empty aux array
            aux[j] = 0.0;
        }
    }
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

    for (int j = 0; j < 3; j++)
    {
        jacobi_x[j] = x_cm[j] / eta[objects_count - 1];
        jacobi_v[j] = v_cm[j] / eta[objects_count - 1];
    }
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

    for (int j = 0; j < 3; j++)
    {
        x[j] =  x_cm[j] / m[0];
        v[j] =  v_cm[j] / m[0];
    }
}

void stumpff_functions(
    real psi,
    real *restrict c2,
    real *restrict c3
)
{
    if (psi > 1e-10)
    {
        real sqrt_psi = sqrt(psi);
        *c2 = (1.0 - cos(sqrt_psi)) / psi;
        *c3 = (sqrt_psi - sin(sqrt_psi)) / (psi * sqrt_psi);
    }
    else if (psi < -1e-10)
    {
        real sqrt_abs_psi = sqrt(-psi);
        *c2 = (1.0 - cosh(sqrt_abs_psi)) / psi;
        *c3 = (sqrt_abs_psi - sinh(sqrt_abs_psi)) / (psi * sqrt_abs_psi);
    }
    else
    {
        *c2 = 0.5;
        *c3 = 1.0 / 6.0;
    }
}

void propagate_kepler(
    int i,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    real gm,
    real dt
)
{
    real temp_vec[3];

    // Energy tolerance: used to distinguish between elliptic, parabolic, and hyperbolic orbits,
    // ideally 0:
    real tol_energy = 0.0;

    real x[3];
    real v[3];
    memcpy(x, &jacobi_x[i * 3], 3 * sizeof(real));
    memcpy(v, &jacobi_v[i * 3], 3 * sizeof(real));
    real x_dot_v = vec_dot(x, v, 3);
    real x_norm = vec_norm(x, 3);
    real v_norm = vec_norm(v, 3);

    // Initial value of the Keplerian energy:
    real xi = v_norm * v_norm * 0.5 - gm / x_norm;

    real semi_major_axis = - gm / (2.0 * xi);
    real alpha = 1.0 / semi_major_axis;
    real sqrt_gm = sqrt(gm);
    real chi_0;

    // Elliptic orbits
    if (alpha > tol_energy + 1e-12)
    {
        chi_0 = sqrt_gm * dt * alpha;
    } 

    // Hyperbolic orbits
    else if (alpha < tol_energy - 1e-12)
    {
        chi_0 = ( 
            sqrt(-semi_major_axis) * log(
                -2.0 * gm * alpha * dt / (
                    x_dot_v + sqrt(-gm * semi_major_axis) * (1.0 - x_norm * alpha)
                )
            )
        ) /* * copysign(1.0, dt) */;
    }

    // Parabolic orbits
    else
    {
        vec_cross(&jacobi_x[i * 3], &jacobi_v[i * 3], temp_vec);
        real p = vec_norm(temp_vec, 3);
        p = p * p / gm;
        real s = 0.5 * atan(1.0 / (3.0 * sqrt(gm / (p * p * p)) * dt));
        real w = atan(pow(tan(s), (1.0 / 3.0)));
        chi_0 = sqrt(p) * 2.0 / tan(2.0 * w); 
    }

    // Solve Kepler's equation
    real psi;
    real c2;
    real c3;
    real r;
    real chi;
    int j = 0;
    for (; j < 500; j++)
    {
        // Compute universal variable
        psi = chi_0 * chi_0 * alpha;

        // Copmute Stumpff functions
        stumpff_functions(psi, &c2, &c3);

        // Propagate radial distance
        r = (
            chi_0 * chi_0 * c2 
            + x_dot_v / sqrt_gm * chi_0 * (1.0 - psi * c3) 
            + x_norm * (1.0 - psi * c2)
        );

        // Auxiliary variable for f and g functions
        chi = chi_0 + (
                sqrt_gm * dt 
                - (chi_0 * chi_0 * chi_0) * c3
                - x_dot_v / sqrt_gm * (chi_0 * chi_0) * c2
                - x_norm * chi_0 * (1.0 - psi * c3)
            ) / r;

        // Check convergence
        if (fabs(chi - chi_0) < TOL_KEPLER)
        {
            break;
        }

        chi_0 = chi;
    }

    if (j == 500)
    {
        printf("Warning: Kepler's equation did not converge, error = %23.15e\n",
            fabs(chi - chi_0)
        );
    }

    // Evaluate f and g functions, together with their derivatives
    real f = 1.0 - (chi * chi) * c2 / x_norm;
    real g = dt - (chi * chi * chi) * c3 / sqrt_gm;

    real df = sqrt_gm / (r * x_norm) * chi * (psi * c3 - 1.0);
    real dg = 1.0 - (chi * chi) * c2 / r;

    // Compute position and velocity vectors
    for (int j = 0; j < 3; j++)
    {
        jacobi_x[i * 3 + j] = f * x[j] + g * v[j];
        jacobi_v[i * 3 + j] = df * x[j] + dg * v[j];
    }
}