/**
 * \file integrator_whfast.c
 * \author Ching Yin Ng
 * \brief Function definitions for WHFast integrators
 * 
 * Function definitions for WHfast integrators. This is
 * a C implementation with modifications based on the reference:
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
#include "math_functions.h"
#include "storing.h"

/**
 * \brief Compute the velocity kick
 * 
 * \param objects_count Number of objects in the system
 * \param jacobi_v Array of Jacobi velocity vectors
 * \param a Array of acceleration vectors
 * \param dt Time step of the system
 */
IN_FILE void whfast_kick(
    const int objects_count,
    real *restrict jacobi_v,
    const real *restrict a,
    const real dt
);

/** 
 * \brief Compute the position drift
 * 
 * \param system Pointer to the gravitational system
 * \param jacobi_x Array of Jacobi position vectors
 * \param jacobi_v Array of Jacobi velocity vectors
 * \param eta Array of cumulative masses
 * \param dt Time step of the system
 * \param kepler_tol Tolerance for solving Kepler's equation
 * \param kepler_max_iter Maximum number of iterations in solving Kepler's equation
 * \param kepler_auto_remove Flag to indicate whether to remove objects
 *                           that failed to converge in Kepler's equation
 * \param kepler_failed_bool_array Array of flags to indicate whether 
 *                                 an object failed to converge in Kepler's equation
 * \param kepler_failed_flag Flag to indicate whether any object failed to converge
 *                           in Kepler's equation
 * \param verbose Verbosity level
 */
IN_FILE int whfast_drift(
    System *restrict system,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    const real *restrict eta,
    const real dt,
    const real kepler_tol,
    const int kepler_max_iter,
    const bool kepler_auto_remove,
    const real kepler_auto_remove_tol,
    bool *restrict kepler_failed_bool_array,
    bool *restrict kepler_failed_flag,
    const int verbose
);

/**
 * \brief Transform Cartesian coordinates to Jacobi coordinates
 * 
 * \param system Pointer to the gravitational system
 * \param jacobi_x Array of Jacobi position vectors to be stored
 * \param jacobi_v Array of Jacobi velocity vectors to be stored
 * \param eta Array of cumulative masses
 */
IN_FILE void cartesian_to_jacobi(
    System *restrict system,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    const real *restrict eta
);

/**
 * \brief Transform Jacobi coordinates to Cartesian coordinates
 * 
 * \param system Pointer to the gravitational system
 * \param jacobi_x Array of Jacobi position vectors
 * \param jacobi_v Array of Jacobi velocity vectors
 * \param eta Array of cumulative masses
 */
IN_FILE void jacobi_to_cartesian(
    System *restrict system,
    const real *restrict jacobi_x,
    const real *restrict jacobi_v,
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
 * 
 * \retval SUCCESS If exit successfully
 * \retval ERROR_WHFAST_STUMPFF_Z_INFINITE If z is infinite
 */
IN_FILE int stumpff_functions(
    real z,
    real *restrict c0,
    real *restrict c1,
    real *restrict c2,
    real *restrict c3
);

/**
 * \brief Direct pairwise acceleration function for WHFast integrator
 * 
 * \details This is a brute-force pairwise calculation
 *          of gravitational acceleration between all objects,
 *          which is O(n^2) complexity.
 * 
 * \param a Array of acceleration vectors to be stored
 * \param system Pointer to the gravitational system
 * \param jacobi_x Array of Jacobi position vectors
 * \param eta Array of cumulative masses
 * \param acceleration_param Pointer to acceleration parameters
 * 
 * \retval SUCCESS If exit successfully
 */
IN_FILE int whfast_acceleration_pairwise(
    real *restrict a,
    const System *system,
    real *restrict jacobi_x,
    const real *restrict eta,
    const AccelerationParam *acceleration_param
);

/**
 * \brief Acceleration function for WHFast integrator,
 *        separating massive and massless objects
 * 
 * \details This function calculates the gravitational acceleration
 *          between massive and massless objects separately.
 *          This is an O(m^2 + mn) complexity calculation,
 *          where m and n are the number of massive and massless 
 *          objects, respectively.
 * 
 * \param a Array of acceleration vectors to be stored
 * \param system Pointer to the gravitational system
 * \param jacobi_x Array of Jacobi position vectors
 * \param eta Array of cumulative masses
 * \param acceleration_param Pointer to acceleration parameters
 * 
 * \retval SUCCESS If exit successfully
 * \retval ERROR_WHFAST_ACCELERATION_MASSLESS_MEMORY_ALLOC If memory allocation failed
 */

IN_FILE int whfast_acceleration_massless(
    real *restrict a,
    const System *system,
    real *restrict jacobi_x,
    const real *restrict eta,
    const AccelerationParam *acceleration_param
);

WIN32DLL_API int whfast(
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

    IN_FILE int (*whfast_acceleration)(
        real *restrict a,
        const System *system,
        real *restrict jacobi_x,
        const real *restrict eta,
        const AccelerationParam *acceleration_param
    );

    if (strcmp(acceleration_param->method, "pairwise") == 0)
    {
        whfast_acceleration = whfast_acceleration_pairwise;
    }
    else if (strcmp(acceleration_param->method, "massless") == 0)
    {
        whfast_acceleration = whfast_acceleration_massless;
    }
    else
    {
        return_code = ERROR_WHFAST_UNKNOWN_ACCELERATION_METHOD;
        goto err_unknown_acc_method;
    }

    int objects_count = system->objects_count;
    real *restrict m = system->m;
    
    const real dt = integrator_param->dt;
    const int64 n_steps = simulation_param->n_steps_;
    const int storing_freq = storing_param->storing_freq;

    const real kepler_tol = integrator_param->whfast_kepler_tol;
    const int kepler_max_iter = integrator_param->whfast_kepler_max_iter;
    const bool kepler_auto_remove = integrator_param->whfast_kepler_auto_remove;
    const real kepler_auto_remove_tol = integrator_param->whfast_kepler_auto_remove_tol;

    /* Allocate memory for calculation */
    real *restrict jacobi_x = calloc(objects_count * 3, sizeof(real));
    real *restrict jacobi_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict temp_jacobi_v = malloc(objects_count * 3 * sizeof(real));
    real *restrict a = malloc(objects_count * 3 * sizeof(real));
    real *restrict eta = malloc(objects_count * sizeof(real));

    if (
        !jacobi_x
        || !jacobi_v
        || !temp_jacobi_v
        || !a
        || !eta
    )
    {
        return_code = ERROR_WHFAST_MEMORY_ALLOC;
        goto err_memory_alloc;
    }

    /* Auto remove objects that failed to converge in Kepler's equation */
    bool *kepler_failed_bool_array = NULL;
    bool kepler_failed_flag = false;
    if (kepler_auto_remove)
    {
        kepler_failed_bool_array = calloc(objects_count, sizeof(bool));

        if (!kepler_failed_bool_array)
        {
            return_code = ERROR_WHFAST_KEPLER_AUTO_REMOVE_MEMORY_ALLOC;
            goto err_kepler_auto_remove_memory;
        }
    }

    /* Initialization */
    eta[0] = m[0];
    for (int i = 1; i < objects_count; i++)
    {
        eta[i] = eta[i - 1] + m[i];
    }
    cartesian_to_jacobi(system, jacobi_x, jacobi_v, eta);
    return_code = whfast_acceleration(a, system, jacobi_x, eta, acceleration_param);
    if (return_code != SUCCESS)
    {
        goto err_acc;
    }
    whfast_kick(objects_count, jacobi_v, a, 0.5 * dt);
    
    /* Main Loop */
    for (int64 count = 1; count <= n_steps; count++)
    {   
        return_code = whfast_drift(
            system,
            jacobi_x,
            jacobi_v,
            eta,
            dt,
            kepler_tol,
            kepler_max_iter,
            kepler_auto_remove,
            kepler_auto_remove_tol,
            kepler_failed_bool_array,
            &kepler_failed_flag,
            settings->verbose
        );
        if (return_code != SUCCESS)
        {
            goto err_drift;
        }
        
        /**
         * Remove objects that failed to converge in Kepler's equation
         * 
         * IMPORTANT:
         * It is important to remove objects right after drift step.
         * Otherwise, one nan value could pollute the data of the 
         * whole system with nan values, especially for large N.
         */
        if (kepler_auto_remove && kepler_failed_flag)
        {
            kepler_failed_flag = false;

            /* Remove object */
            int kepler_remove_count = 0;

            // The first object is the central object, which is not 
            // calculated in drift step
            for (int i = 1; i < objects_count; i++)
            {
                if (kepler_failed_bool_array[i])
                {
                    kepler_failed_bool_array[i] = false;
                    kepler_remove_count++;
                    if (settings->verbose > 0)
                    {
                        fprintf(stderr, "kepler_auto_remove: Object %d with mass %f removed\n", i, m[i]);
                    }
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
            system->objects_count = objects_count;
            if (settings->verbose > 0)
            {
                fprintf(stderr, "kepler_auto_remove: %d objects removed in total. \
                                 Remaining objects: %d\n", kepler_remove_count,
                                 objects_count);
            }
        }
        jacobi_to_cartesian(system, jacobi_x, jacobi_v, eta);
        return_code = whfast_acceleration(a, system, jacobi_x, eta, acceleration_param);
        if (return_code != SUCCESS)
        {
            goto err_acc;
        }
        whfast_kick(objects_count, jacobi_v, a, dt);

        *(simulation_status->t) = count * dt;

        /* Store solution */
        if (count % storing_freq == 0)
        {
            // Get v_1 from v_1+1/2
            memcpy(temp_jacobi_v, jacobi_v, objects_count * 3 * sizeof(real));
            whfast_kick(objects_count, temp_jacobi_v, a, -0.5 * dt);
            jacobi_to_cartesian(system, jacobi_x, jacobi_v, eta);
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

    /* Free memory */
    free(jacobi_x);
    free(jacobi_v);
    free(temp_jacobi_v);
    free(a);
    free(eta);
    free(kepler_failed_bool_array);

    return SUCCESS;

err_user_interrupt:
err_store_solution:
err_acc:
err_drift:
err_kepler_auto_remove_memory:
    free(kepler_failed_bool_array);
err_memory_alloc:
    free(eta);
    free(a);
    free(temp_jacobi_v);
    free(jacobi_v);
    free(jacobi_x);
err_unknown_acc_method:
    return return_code;
}

IN_FILE void whfast_kick(
    const int objects_count,
    real *restrict jacobi_v,
    const real *restrict a,
    const real dt
)
{
    for (int i = 0; i < objects_count; i++)
    {
        jacobi_v[i * 3 + 0] += a[i * 3 + 0] * dt;
        jacobi_v[i * 3 + 1] += a[i * 3 + 1] * dt;
        jacobi_v[i * 3 + 2] += a[i * 3 + 2] * dt;
    }
}

IN_FILE int whfast_drift(
    System *restrict system,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    const real *restrict eta,
    const real dt,
    const real kepler_tol,
    const int kepler_max_iter,
    const bool kepler_auto_remove,
    const real kepler_auto_remove_tol,
    bool *restrict kepler_failed_bool_array,
    bool *restrict kepler_failed_flag,
    const int verbose
)
{
    int return_code;

    const int objects_count = system->objects_count;
    const real *restrict m = system->m;
    const real G = system->G;
    for (int i = 1; i < objects_count; i++)
    {
        real gm = G * m[0] * eta[i] / eta[i - 1];
        real x[3];
        real v[3];
        memcpy(x, &jacobi_x[i * 3], 3 * sizeof(real));
        memcpy(v, &jacobi_v[i * 3], 3 * sizeof(real));

        real x_norm = vec_norm_3d(x);
        real v_norm = vec_norm_3d(v);

        // Radial velocity
        real radial_v = vec_dot_3d(x, v) / x_norm; 

        real alpha = 2.0 * gm / x_norm - (v_norm * v_norm);

        /* Solve Kepler's equation with Newton-Raphson method */
        
        // Initial guess
        real s = dt / x_norm;

        // Solve Kepler's equation
        real c0 = 0.0;
        real c1 = 0.0;
        real c2 = 0.0;
        real c3 = 0.0;
        bool is_converged = false;

        for (int j = 0; j < kepler_max_iter; j++)
        {
            // Compute Stumpff functions
            return_code = stumpff_functions(alpha * (s * s), &c0, &c1, &c2, &c3);
            if (return_code != SUCCESS)
            {
                goto err_stumpff;
            }

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
                is_converged = true;
                break;
            }
        }

        // The raidal distance is equal to the derivative of F
        // real r = dF
        real r = x_norm * c0 + x_norm * radial_v * s * c1 + gm * (s * s) * c2;

        if (!is_converged)
        {
            real error = (
                x_norm * s * c1
                + x_norm * radial_v * (s * s) * c2
                + gm * (s * s * s) * c3
                - dt
            ) / r;
            
            if (verbose >= 2)
            {
                fprintf(stderr, "Warning: Kepler's equation did not converge. \
                                Object index: %d, error = %23.15g\n", i, error);
            }

            if (kepler_auto_remove && ((fabs(error) > kepler_auto_remove_tol) || isnan(error)))
            {
                kepler_failed_bool_array[i] = true;
                *kepler_failed_flag = true;
            }
        }

        /* Evaluate f and g functions, together with their derivatives */
        real f = 1.0 - gm * (s * s) * c2 / x_norm;
        real g = dt - gm * (s * s * s) * c3;

        real df = -gm * s * c1 / (r * x_norm);
        real dg = 1.0 - gm * (s * s) * c2 / r; 

        /* Compute position and velocity vectors */
        for (int j = 0; j < 3; j++)
        {
            jacobi_x[i * 3 + j] = f * x[j] + g * v[j];
            jacobi_v[i * 3 + j] = df * x[j] + dg * v[j];
        }
    }

    return SUCCESS;

err_stumpff:
    return return_code;
}

IN_FILE void cartesian_to_jacobi(
    System *restrict system,
    real *restrict jacobi_x,
    real *restrict jacobi_v,
    const real *restrict eta
)
{
    real x_cm[3];
    real v_cm[3];
    const int objects_count = system->objects_count;
    real *restrict x = system->x;
    real *restrict v = system->v;
    const real *restrict m = system->m;

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

IN_FILE void jacobi_to_cartesian(
    System *restrict system,
    const real *restrict jacobi_x,
    const real *restrict jacobi_v,
    const real *restrict eta
)
{
    real x_cm[3];
    real v_cm[3];
    const int objects_count = system->objects_count;
    real *restrict x = system->x;
    real *restrict v = system->v;
    const real *restrict m = system->m;

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

IN_FILE int stumpff_functions(
    real z,
    real *restrict c0,
    real *restrict c1,
    real *restrict c2,
    real *restrict c3
)
{
    int return_code;

    if (isinf(z))
    {
        return_code = ERROR_WHFAST_STUMPFF_Z_INFINITE;
        goto err_z_inf;
    }
    else if (isnan(z))
    {
        return_code = ERROR_WHFAST_STUMPFF_Z_NAN;
        goto err_z_nan;
    }

    /* Reduce the argument */
    int n = 0;
    while (fabs(z) > 0.1)
    {
        z /= 4.0;
        n++;
    }

    /* Compute stumpff functions */
    real temp_c3 = (
        1.0 - z / 20.0 * (1.0 - z / 42.0 * (1.0 - z / 72.0 * (1.0 - z / 110.0 \
        * (1.0 - z / 156.0 * (1.0 - z / 210.0)))))
    ) / 6.0;
    real temp_c2 = (
        1.0 - z / 12.0 * (1.0 - z / 30.0 * (1.0 - z / 56.0 * (1.0 - z / 90.0 \
        * (1.0 - z / 132.0 * (1.0 - z / 182.0)))))
    ) / 2.0;
    real temp_c1 = 1.0 - z * temp_c3;
    real temp_c0 = 1.0 - z * temp_c2;

    /* Half-angle formulae to recover the actual argument */
    while (n > 0)
    {
        temp_c3 = (temp_c2 + temp_c0 * temp_c3) / 4.0;
        temp_c2 = (temp_c1 * temp_c1) / 2.0;
        temp_c1 = temp_c0 * temp_c1;
        temp_c0 = (2.0 * temp_c0 * temp_c0) - 1.0;
        n--;
    }
    *c3 = temp_c3;
    *c2 = temp_c2;
    *c1 = temp_c1;
    *c0 = temp_c0;

    return SUCCESS;

err_z_nan:
err_z_inf:
    return return_code;
}

IN_FILE int whfast_acceleration_pairwise(
    real *restrict a,
    const System *system,
    real *restrict jacobi_x,
    const real *restrict eta,
    const AccelerationParam *acceleration_param
)
{
    const int objects_count = system->objects_count;
    const real *restrict x = system->x;
    const real *restrict m = system->m;
    const real G = system->G;

    const real softening_length = acceleration_param->softening_length;

    real aux[3];
    real temp_vec[3];
    real temp_vec_norm;
    real temp_vec_norm_cube;
    real temp_jacobi_norm;
    real temp_jacobi_norm_cube;
    real softening_length_cube = softening_length * softening_length * softening_length;
    for (int i = 1; i < objects_count; i++)
    {
        // Calculate x_0i
        temp_vec[0] = x[i * 3 + 0] - x[0];
        temp_vec[1] = x[i * 3 + 1] - x[1];
        temp_vec[2] = x[i * 3 + 2] - x[2];

        temp_vec_norm = vec_norm_3d(temp_vec);
        temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;
        temp_jacobi_norm = vec_norm_3d(&jacobi_x[i * 3]);
        temp_jacobi_norm_cube = (temp_jacobi_norm * temp_jacobi_norm * temp_jacobi_norm) + softening_length_cube;
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

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

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

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

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

                temp_vec_norm = vec_norm_3d(temp_vec);
                temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

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

    return SUCCESS;
}

IN_FILE int whfast_acceleration_massless(
    real *restrict a,
    const System *system,
    real *restrict jacobi_x,
    const real *restrict eta,
    const AccelerationParam *acceleration_param
)
{
    int return_code;

    const int objects_count = system->objects_count;
    const real *restrict x = system->x;
    const real *restrict m = system->m;
    const real G = system->G;

    const real softening_length = acceleration_param->softening_length;

    real aux[3];
    real temp_vec[3];
    real temp_vec_norm;
    real temp_vec_norm_cube;
    real temp_jacobi_norm;
    real temp_jacobi_norm_cube;
    real softening_length_cube = softening_length * softening_length * softening_length;

    /* Find the numbers of massive and massless objects */
    int massive_objects_count = 0;
    int massless_objects_count = 0;
    for (int i = 0; i < objects_count; i++)
    {
        if (m[i] != 0.0)
        {
            massive_objects_count++;
        }
        else
        {
            massless_objects_count++;
        }
    }

    /* Find the indices of massive and massless objects */
    int *restrict massive_indices = malloc(massive_objects_count * sizeof(int));
    int *restrict massless_indices = malloc(massless_objects_count * sizeof(int));
    massive_objects_count = 0;
    massless_objects_count = 0;

    if (!massive_indices || !massless_indices)
    {
        return_code = ERROR_WHFAST_ACC_MASSLESS_MEMORY_ALLOC;
        goto err_memory;
    }

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

    /* Acceleration calculation for massive objects */
    for (int i = 1; i < massive_objects_count; i++)
    {
        int idx_i = massive_indices[i];

        // Calculate x_0i
        temp_vec[0] = x[idx_i * 3 + 0] - x[0];
        temp_vec[1] = x[idx_i * 3 + 1] - x[1];
        temp_vec[2] = x[idx_i * 3 + 2] - x[2];

        temp_vec_norm = vec_norm_3d(temp_vec);
        temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;
        temp_jacobi_norm = vec_norm_3d(&jacobi_x[idx_i * 3]);
        temp_jacobi_norm_cube = (temp_jacobi_norm * temp_jacobi_norm * temp_jacobi_norm) + softening_length_cube;
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
            int idx_j = massive_indices[j];

            // Calculate x_ji
            temp_vec[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            temp_vec[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            temp_vec[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

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
            int idx_j = massive_indices[j];

            // Calculate x_ij
            temp_vec[0] = x[idx_j * 3 + 0] - x[idx_i * 3 + 0];
            temp_vec[1] = x[idx_j * 3 + 1] - x[idx_i * 3 + 1];
            temp_vec[2] = x[idx_j * 3 + 2] - x[idx_i * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

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
            int idx_j = massive_indices[j];

            for (int k = i + 1; k < massive_objects_count; k++)
            {
                int idx_k = massive_indices[k];

                // Calculate x_jk
                temp_vec[0] = x[idx_k * 3 + 0] - x[idx_j * 3 + 0];
                temp_vec[1] = x[idx_k * 3 + 1] - x[idx_j * 3 + 1];
                temp_vec[2] = x[idx_k * 3 + 2] - x[idx_j * 3 + 2];

                temp_vec_norm = vec_norm_3d(temp_vec);
                temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

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

    /* Acceleration calculation for massless objects */
    for (int i = 0; i < massless_objects_count; i++)
    {
        int idx_i = massless_indices[i];
        if (idx_i == 0)
        {
            continue;
        }
        
        // Calculate x_0i
        temp_vec[0] = x[idx_i * 3 + 0] - x[0];
        temp_vec[1] = x[idx_i * 3 + 1] - x[1];
        temp_vec[2] = x[idx_i * 3 + 2] - x[2];

        temp_vec_norm = vec_norm_3d(temp_vec);
        temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;
        temp_jacobi_norm = vec_norm_3d(&jacobi_x[idx_i * 3]);
        temp_jacobi_norm_cube = (temp_jacobi_norm * temp_jacobi_norm * temp_jacobi_norm) + softening_length_cube;
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
            int idx_j = massive_indices[j];
            if (idx_j >= idx_i)
            {
                break;
            }

            // Calculate x_ji
            temp_vec[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            temp_vec[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            temp_vec[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

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
            int idx_j = massive_indices[j];
            if (idx_j <= idx_i)
            {
                continue;
            }

            // Calculate x_ij
            temp_vec[0] = x[idx_j * 3 + 0] - x[idx_i * 3 + 0];
            temp_vec[1] = x[idx_j * 3 + 1] - x[idx_i * 3 + 1];
            temp_vec[2] = x[idx_j * 3 + 2] - x[idx_i * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

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
            int idx_j = massive_indices[j];
            if (idx_j >= idx_i)
            {
                break;
            }

            for (int k = j + 1; k < massive_objects_count; k++)
            {
                int idx_k = massive_indices[k];
                if (idx_k <= idx_i)
                {
                    continue;
                }

                // Calculate x_jk
                temp_vec[0] = x[idx_k * 3 + 0] - x[idx_j * 3 + 0];
                temp_vec[1] = x[idx_k * 3 + 1] - x[idx_j * 3 + 1];
                temp_vec[2] = x[idx_k * 3 + 2] - x[idx_j * 3 + 2];

                temp_vec_norm = vec_norm_3d(temp_vec);
                temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

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

    return SUCCESS;

err_memory:
    free(massive_indices);
    free(massless_indices);
    return return_code;
}
