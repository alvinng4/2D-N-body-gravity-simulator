#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
// #include <stdio.h> // For testing

#ifdef WIN32DLL_EXPORTS
    #define WIN32DLL_API __declspec(dllexport)
#else
    #define WIN32DLL_API 
#endif

typedef int64_t int64;
typedef double real;

real abs_max_vec(const real *restrict vec, int vec_length);
real abs_max_vec_array(real (*restrict arr)[3], int objects_count);
real vec_norm(const real *restrict vec, int vec_length);
void compute_energy(
    int objects_count, 
    int npts,
    int *restrict count, 
    real *restrict energy, 
    const real (*restrict sol_state)[objects_count * 6], 
    const real *restrict m, 
    real G
);
void acceleration(int objects_count, real (*restrict x)[3], real (*restrict a)[3], const real *restrict m, real G);
void euler(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    int64 npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    int64 *restrict count,
    int *restrict store_count,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3]
);
void euler_cromer(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    int64 npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    int64 *restrict count,
    int *restrict store_count,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3]
);
void rk4(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    int64 npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    int64 *restrict count,
    int *restrict store_count,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3]
);
void leapfrog(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    int64 npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    int64 *restrict count,
    int *restrict store_count,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3]
);
int rk_embedded(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real *restrict t, 
    real *restrict dt, 
    real tf, 
    int store_every_n,
    int *restrict store_count,
    int *restrict count, 
    const real *restrict m, 
    real G, 
    int power,
    int power_test,
    int len_coeff,
    const real (*restrict coeff)[len_coeff],
    int len_weights,
    const real *restrict weights,
    const real *restrict weights_test,
    real abs_tolerance,
    real rel_tolerance,
    real (*restrict sol_state)[6 * objects_count],
    int len_sol_time,
    real *restrict sol_time,
    real *restrict sol_dt,
    real (*restrict x_err_comp_sum)[3],
    real (*restrict v_err_comp_sum)[3]
);
int ias15(
    int objects_count, 
    int dim_nodes,
    const real *restrict nodes,
    const real *restrict aux_c, 
    const real *restrict aux_r,
    real *restrict aux_b0,
    real *restrict aux_b,
    real *restrict aux_g,
    real *restrict aux_e,
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real (*restrict a)[3],
    const real *restrict m, 
    real G,     
    real *restrict t, 
    real *restrict dt, 
    real tf, 
    int store_every_n,
    int *restrict store_count,
    int *restrict count, 
    real tolerance,
    real tolerance_pc,
    real (*restrict sol_state)[6 * objects_count],
    int len_sol_time,
    real *restrict sol_time,
    real *restrict sol_dt,
    real safety_fac,
    real exponent,
    int *restrict ias15_refine_flag,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3]
);
void ias15_step(
    int objects_count,
    int dim_nodes,
    int dim_nodes_minus_1,
    int dim_nodes_minus_2,
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    const real *restrict m,
    real G,
    real *restrict t,
    real *restrict dt,
    real tf,
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
    real *restrict delta_aux_b,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3],
    real (*restrict temp_x_err_comp_sum)[3], 
    real (*restrict temp_v_err_comp_sum)[3]
);
void ias15_approx_pos_aux(
    int objects_count,
    real (*restrict x)[3],
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real node,
    real *restrict aux_b,
    real dt,
    real (*restrict x_err_comp_sum)[3]
);
void ias15_approx_vel_aux(
    int objects_count,
    real (*restrict v)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real node,
    real *restrict aux_b,
    real dt,
    real (*restrict v_err_comp_sum)[3]
);
void ias15_approx_pos_step(
    int objects_count,
    real (*restrict x)[3],
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real *restrict aux_b,
    real dt,
    real (*restrict temp_x_err_comp_sum)[3]
);
void ias15_approx_vel_step(
    int objects_count,
    real (*restrict v)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real *restrict aux_b,
    real dt,
    real (*restrict temp_v_err_comp_sum)[3]
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

// Find the max absolute value in a 1D array
WIN32DLL_API real abs_max_vec(const real *restrict vec, int vec_length)
{
    real max = fabs(vec[0]);
    for (int i = 1; i < vec_length; i++)
    {
        max = fmax(max, fabs(vec[i]));
    }

    return max;
}

// Find the max absolute value in a array with 3D vectors
WIN32DLL_API real abs_max_vec_array(real (*restrict arr)[3], int objects_count)
{
    real max = 0;
    for (int i = 0; i < objects_count; i++)
    {
        max = fmax(max, fabs(arr[i][0]));
        max = fmax(max, fabs(arr[i][1]));
        max = fmax(max, fabs(arr[i][2]));
    }

    return max;
}

// Find the norm of 1D vector
WIN32DLL_API real vec_norm(const real *restrict vec, int vec_length)
{   
    real sum = 0.0;
    if (vec_length == 3) 
    {
        sum = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    }
    else
    {
        for (int i = 0; i < vec_length; i++) 
        {
            sum += vec[i] * vec[i];
        }
    }
    return sqrt(sum);
}

WIN32DLL_API void compute_energy(
    int objects_count, 
    int npts,
    int *restrict count, 
    real *restrict energy, 
    const real (*restrict sol_state)[objects_count * 6], 
    const real *restrict m, 
    real G
)
{
    // Round down current progress percentage as int
    int progress_percentage = (*count / npts * 100);

    real temp_vec[3];

    while (1)
    {   
        for (int i = 0; i < objects_count; i++)
        {
            // KE
            energy[*count] += (
                0.5 * m[i] 
                * pow(vec_norm(&sol_state[*count][(objects_count + i) * 3], 3), 2)
            );

            // PE
            for (int j = i + 1; j < objects_count; j++)
            {
                temp_vec[0] = (
                    sol_state[*count][i * 3 + 0] 
                    - sol_state[*count][j * 3 + 0]
                );
                temp_vec[1] = (
                    sol_state[*count][i * 3 + 1] 
                    - sol_state[*count][j * 3 + 1]
                );
                temp_vec[2] = (
                    sol_state[*count][i * 3 + 2] 
                    - sol_state[*count][j * 3 + 2]
                );
                energy[*count] -= (
                    G * m[i] * m[j]
                    / vec_norm(temp_vec, 3)
                );
            }
        }
        *count += 1;

        if (*count >= npts)
        {
            break;
        }

        // Exit to update progress bar
        if ((*count / npts * 100) > progress_percentage)
        {   
            break;
        }
    }
}

WIN32DLL_API void acceleration(int objects_count, real (*restrict x)[3], real (*restrict a)[3], const real *restrict m, real G)
{   
    real R_norm;
    real temp_value;
    real temp_vec[3];
    real R[3];

    // Empty the input array
    for (int i = 0; i < objects_count; i++)
    {
        a[i][0] = 0.0;
        a[i][1] = 0.0;
        a[i][2] = 0.0;
    }

    for(int i = 0; i < objects_count; i++)
    {
        for(int j = i + 1; j < objects_count; j++)
        {
            // Calculate \vec{R} and its norm
            R[0] = x[i][0] - x[j][0];
            R[1] = x[i][1] - x[j][1];
            R[2] = x[i][2] - x[j][2];
            R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[i][0] -= temp_vec[0] * m[j];
            a[i][1] -= temp_vec[1] * m[j];
            a[i][2] -= temp_vec[2] * m[j];
            a[j][0] += temp_vec[0] * m[i];
            a[j][1] += temp_vec[1] * m[i];
            a[j][2] += temp_vec[2] * m[i];
        }
    }
}

WIN32DLL_API void euler(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    int64 npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    int64 *restrict count,
    int *restrict store_count,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3]
)
{   
    real (*temp_x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    // Round down current progress percentage as int
    int progress_percentage = *count * 100 / npts;

    // Main Loop
    while ((*count + 1) <= npts)
    {   
        acceleration(objects_count, x, a, m, G);

        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            x_err_comp_sum[j][0] += v[j][0] * dt;
            x_err_comp_sum[j][1] += v[j][1] * dt;
            x_err_comp_sum[j][2] += v[j][2] * dt;
            v_err_comp_sum[j][0] += a[j][0] * dt;
            v_err_comp_sum[j][1] += a[j][1] * dt;
            v_err_comp_sum[j][2] += a[j][2] * dt;

            x[j][0] = temp_x[j][0] + x_err_comp_sum[j][0];
            x[j][1] = temp_x[j][1] + x_err_comp_sum[j][1];
            x[j][2] = temp_x[j][2] + x_err_comp_sum[j][2];
            v[j][0] = temp_v[j][0] + v_err_comp_sum[j][0];
            v[j][1] = temp_v[j][1] + v_err_comp_sum[j][1];
            v[j][2] = temp_v[j][2] + v_err_comp_sum[j][2];

            x_err_comp_sum[j][0] += temp_x[j][0] - x[j][0];
            x_err_comp_sum[j][1] += temp_x[j][1] - x[j][1];
            x_err_comp_sum[j][2] += temp_x[j][2] - x[j][2];
            v_err_comp_sum[j][0] += temp_v[j][0] - v[j][0];
            v_err_comp_sum[j][1] += temp_v[j][1] - v[j][1];
            v_err_comp_sum[j][2] += temp_v[j][2] - v[j][2];

            // Store solution
            if ((*count + 1) % store_every_n == 0)
            {
                sol_state[*store_count + 1][j * 3] = x[j][0];
                sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }

            if ((*count + 2) == npts)
            {
                sol_state[store_npts - 1][j * 3] = x[j][0];
                sol_state[store_npts - 1][j * 3 + 1] = x[j][1];
                sol_state[store_npts - 1][j * 3 + 2] = x[j][2];
                sol_state[store_npts - 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }
        }    
        if ((*count + 1) % store_every_n == 0)
        {
            *store_count += 1;
        }

        *count += 1;

        // Exit to update progress bar
        if ((*count * 100 / npts) > progress_percentage)
        {   
            free(a);
            free(temp_x);
            free(temp_v);
            return;
        }
    }
}

WIN32DLL_API void euler_cromer(
    int objects_count,
    real (*restrict x)[3],
    real (*restrict v)[3],
    real dt,
    int64 npts,
    const real *restrict m,
    real G,
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    int64 *restrict count,
    int *restrict store_count,
    real (*restrict x_err_comp_sum)[3],
    real (*restrict v_err_comp_sum)[3]
)
{   
    real (*temp_x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    // Round down current progress percentage as int
    int progress_percentage = *count * 100 / npts;

    // Main Loop
    while ((*count + 1) <= npts)
    {   
        acceleration(objects_count, x, a, m, G);

        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation of v
            v_err_comp_sum[j][0] += a[j][0] * dt;
            v_err_comp_sum[j][1] += a[j][1] * dt;
            v_err_comp_sum[j][2] += a[j][2] * dt;

            v[j][0] = temp_v[j][0] + v_err_comp_sum[j][0];
            v[j][1] = temp_v[j][1] + v_err_comp_sum[j][1];
            v[j][2] = temp_v[j][2] + v_err_comp_sum[j][2];

            v_err_comp_sum[j][0] += temp_v[j][0] - v[j][0];
            v_err_comp_sum[j][1] += temp_v[j][1] - v[j][1];
            v_err_comp_sum[j][2] += temp_v[j][2] - v[j][2];

            // Calculation of x
            x_err_comp_sum[j][0] += v[j][0] * dt;
            x_err_comp_sum[j][1] += v[j][1] * dt;
            x_err_comp_sum[j][2] += v[j][2] * dt;

            x[j][0] = temp_x[j][0] + x_err_comp_sum[j][0];
            x[j][1] = temp_x[j][1] + x_err_comp_sum[j][1];
            x[j][2] = temp_x[j][2] + x_err_comp_sum[j][2];

            x_err_comp_sum[j][0] += temp_x[j][0] - x[j][0];
            x_err_comp_sum[j][1] += temp_x[j][1] - x[j][1];
            x_err_comp_sum[j][2] += temp_x[j][2] - x[j][2];

            // Store solution
            if ((*count + 1) % store_every_n == 0)
            {
                sol_state[*store_count + 1][j * 3] = x[j][0];
                sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }

            if ((*count + 2) == npts)
            {
                sol_state[store_npts - 1][j * 3] = x[j][0];
                sol_state[store_npts - 1][j * 3 + 1] = x[j][1];
                sol_state[store_npts - 1][j * 3 + 2] = x[j][2];
                sol_state[store_npts - 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }
        }    
        if ((*count + 1) % store_every_n == 0)
        {
            *store_count += 1;
        }

        *count += 1;

        // Exit to update progress bar
        if ((*count * 100 / npts) > progress_percentage)
        {   
            free(a);
            free(temp_x);
            free(temp_v);
            return;
        }
    }
}

WIN32DLL_API void rk4(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    int64 npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    int64 *restrict count,
    int *restrict store_count,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3]
)
{
    real (*temp_x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_v)[3] = malloc(objects_count * 3 * sizeof(real));

    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    real (*vk1)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*vk2)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*vk3)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*vk4)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*xk1)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*xk2)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*xk3)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*xk4)[3] = malloc(objects_count * 3 * sizeof(real));

    // Round down current progress percentage as int
    int progress_percentage = *count * 100 / npts;

    // Main Loop
    while ((*count + 1) <= npts)
    {   
        acceleration(objects_count, x, a, m, G);
        memcpy(vk1, a, objects_count * 3 * sizeof(real));
        memcpy(xk1, v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            temp_x[j][0] = x[j][0] + 0.5 * xk1[j][0] * dt;
            temp_x[j][1] = x[j][1] + 0.5 * xk1[j][1] * dt;
            temp_x[j][2] = x[j][2] + 0.5 * xk1[j][2] * dt;
            temp_v[j][0] = v[j][0] + 0.5 * vk1[j][0] * dt;
            temp_v[j][1] = v[j][1] + 0.5 * vk1[j][1] * dt;
            temp_v[j][2] = v[j][2] + 0.5 * vk1[j][2] * dt;
        }
        acceleration(objects_count, temp_x, a, m, G);
        memcpy(vk2, a, objects_count * 3 * sizeof(real));
        memcpy(xk2, temp_v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            temp_x[j][0] = x[j][0] + 0.5 * xk2[j][0] * dt;
            temp_x[j][1] = x[j][1] + 0.5 * xk2[j][1] * dt;
            temp_x[j][2] = x[j][2] + 0.5 * xk2[j][2] * dt;
            temp_v[j][0] = v[j][0] + 0.5 * vk2[j][0] * dt;
            temp_v[j][1] = v[j][1] + 0.5 * vk2[j][1] * dt;
            temp_v[j][2] = v[j][2] + 0.5 * vk2[j][2] * dt;
        }
        acceleration(objects_count, temp_x, a, m, G);
        memcpy(vk3, a, objects_count * 3 * sizeof(real));
        memcpy(xk3, temp_v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            temp_x[j][0] = x[j][0] + xk3[j][0] * dt;
            temp_x[j][1] = x[j][1] + xk3[j][1] * dt;
            temp_x[j][2] = x[j][2] + xk3[j][2] * dt;
            temp_v[j][0] = v[j][0] + vk3[j][0] * dt;
            temp_v[j][1] = v[j][1] + vk3[j][1] * dt;
            temp_v[j][2] = v[j][2] + vk3[j][2] * dt;
        }
        acceleration(objects_count, temp_x, a, m, G);
        memcpy(vk4, a, objects_count * 3 * sizeof(real));
        memcpy(xk4, temp_v, objects_count * 3 * sizeof(real));

        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        memcpy(temp_x, x, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            v_err_comp_sum[j][0] += (vk1[j][0] + 2 * vk2[j][0] + 2 * vk3[j][0] + vk4[j][0]) * dt / 6.0;
            v_err_comp_sum[j][1] += (vk1[j][1] + 2 * vk2[j][1] + 2 * vk3[j][1] + vk4[j][1]) * dt / 6.0;
            v_err_comp_sum[j][2] += (vk1[j][2] + 2 * vk2[j][2] + 2 * vk3[j][2] + vk4[j][2]) * dt / 6.0;
            x_err_comp_sum[j][0] += (xk1[j][0] + 2 * xk2[j][0] + 2 * xk3[j][0] + xk4[j][0]) * dt / 6.0;
            x_err_comp_sum[j][1] += (xk1[j][1] + 2 * xk2[j][1] + 2 * xk3[j][1] + xk4[j][1]) * dt / 6.0;
            x_err_comp_sum[j][2] += (xk1[j][2] + 2 * xk2[j][2] + 2 * xk3[j][2] + xk4[j][2]) * dt / 6.0;

            v[j][0] = temp_v[j][0] + v_err_comp_sum[j][0];
            v[j][1] = temp_v[j][1] + v_err_comp_sum[j][1];
            v[j][2] = temp_v[j][2] + v_err_comp_sum[j][2];
            x[j][0] = temp_x[j][0] + x_err_comp_sum[j][0];
            x[j][1] = temp_x[j][1] + x_err_comp_sum[j][1];
            x[j][2] = temp_x[j][2] + x_err_comp_sum[j][2];

            v_err_comp_sum[j][0] += temp_v[j][0] - v[j][0];
            v_err_comp_sum[j][1] += temp_v[j][1] - v[j][1];
            v_err_comp_sum[j][2] += temp_v[j][2] - v[j][2];
            x_err_comp_sum[j][0] += temp_x[j][0] - x[j][0];
            x_err_comp_sum[j][1] += temp_x[j][1] - x[j][1];
            x_err_comp_sum[j][2] += temp_x[j][2] - x[j][2];

            // Store solution
            if ((*count + 1) % store_every_n == 0)
            {
                sol_state[*store_count + 1][j * 3] = x[j][0];
                sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }

            if ((*count + 2) == npts)
            {
                sol_state[store_npts - 1][j * 3] = x[j][0];
                sol_state[store_npts - 1][j * 3 + 1] = x[j][1];
                sol_state[store_npts - 1][j * 3 + 2] = x[j][2];
                sol_state[store_npts - 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }
        }    
        if ((*count + 1) % store_every_n == 0)
        {
            *store_count += 1;
        }

        *count += 1;

        // Exit to update progress bar
        if ((*count * 100 / npts) > progress_percentage)
        {   
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
            return;
        }
    }
}

WIN32DLL_API void leapfrog(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    int64 npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    int64 *restrict count,
    int *restrict store_count,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3]
)
{   
    real (*temp_x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a_0)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a_1)[3] = malloc(objects_count * 3 * sizeof(real));

    int is_initialize = 1;

    // Round down current progress percentage as int
    int progress_percentage = *count * 100 / npts;

    // Main Loop
    while ((*count + 1) <= npts)
    {       
        if (is_initialize == 1)
        {
            acceleration(objects_count, x, a_0, m, G);
            is_initialize = 0;
        }
        else 
        {
            // Use a_1 from last iteration as a_0
            memcpy(a_0, a_1, objects_count * 3 * sizeof(real));
        }

        // Calculation of x
        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        for (int j = 0; j < objects_count; j++)
        {
            x_err_comp_sum[j][0] += v[j][0] * dt + 0.5 * a_0[j][0] * dt * dt;
            x_err_comp_sum[j][1] += v[j][1] * dt + 0.5 * a_0[j][1] * dt * dt;
            x_err_comp_sum[j][2] += v[j][2] * dt + 0.5 * a_0[j][2] * dt * dt;

            x[j][0] = temp_x[j][0] + x_err_comp_sum[j][0];
            x[j][1] = temp_x[j][1] + x_err_comp_sum[j][1];
            x[j][2] = temp_x[j][2] + x_err_comp_sum[j][2];

            x_err_comp_sum[j][0] += temp_x[j][0] - x[j][0];
            x_err_comp_sum[j][1] += temp_x[j][1] - x[j][1];
            x_err_comp_sum[j][2] += temp_x[j][2] - x[j][2];
        }    

        acceleration(objects_count, x, a_1, m, G);
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation of v
            v_err_comp_sum[j][0] += 0.5 * (a_0[j][0] + a_1[j][0]) * dt;
            v_err_comp_sum[j][1] += 0.5 * (a_0[j][1] + a_1[j][1]) * dt;
            v_err_comp_sum[j][2] += 0.5 * (a_0[j][2] + a_1[j][2]) * dt;

            v[j][0] = temp_v[j][0] + v_err_comp_sum[j][0];
            v[j][1] = temp_v[j][1] + v_err_comp_sum[j][1];
            v[j][2] = temp_v[j][2] + v_err_comp_sum[j][2];

            v_err_comp_sum[j][0] += temp_v[j][0] - v[j][0];
            v_err_comp_sum[j][1] += temp_v[j][1] - v[j][1];
            v_err_comp_sum[j][2] += temp_v[j][2] - v[j][2];

            // Store solution
            if ((*count + 1) % store_every_n == 0)
            {
                sol_state[*store_count + 1][j * 3] = x[j][0];
                sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }

            if ((*count + 2) == npts)
            {
                sol_state[store_npts - 1][j * 3] = x[j][0];
                sol_state[store_npts - 1][j * 3 + 1] = x[j][1];
                sol_state[store_npts - 1][j * 3 + 2] = x[j][2];
                sol_state[store_npts - 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }
        }    
        if ((*count + 1) % store_every_n == 0)
        {
            *store_count += 1;
        }

        *count += 1;

        // Exit to update progress bar
        if ((*count * 100 / npts) > progress_percentage)
        {   
            free(a_0);
            free(a_1);
            free(temp_x);
            free(temp_v);
            return;
        }
    }
}

WIN32DLL_API int rk_embedded(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real *restrict t, 
    real *restrict dt, 
    real tf, 
    int store_every_n,
    int *restrict store_count,
    int *restrict count, 
    const real *restrict m, 
    real G, 
    int power,
    int power_test,
    int len_coeff,
    const real (*restrict coeff)[len_coeff],
    int len_weights,
    const real *restrict weights,
    const real *restrict weights_test,
    real abs_tolerance,
    real rel_tolerance,
    real (*restrict sol_state)[6 * objects_count],
    int len_sol_time,
    real *restrict sol_time,
    real *restrict sol_dt,
    real (*restrict x_err_comp_sum)[3],
    real (*restrict v_err_comp_sum)[3]
)
{
    // Initialization
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
    real (*temp_x_err_comp_sum)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_v_err_comp_sum)[3] = malloc(objects_count * 3 * sizeof(real));

    // Round down current progress percentage as int
    int progress_percentage = (int) (*t / tf * 100.0);

    // Main Loop
    while (1)
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
                temp_v[k][0] = v[k][0] + *dt * temp_v[k][0] + v_err_comp_sum[k][0];
                temp_v[k][1] = v[k][1] + *dt * temp_v[k][1] + v_err_comp_sum[k][1];
                temp_v[k][2] = v[k][2] + *dt * temp_v[k][2] + v_err_comp_sum[k][2];
                temp_x[k][0] = x[k][0] + *dt * temp_x[k][0] + x_err_comp_sum[k][0];
                temp_x[k][1] = x[k][1] + *dt * temp_x[k][1] + x_err_comp_sum[k][1];
                temp_x[k][2] = x[k][2] + *dt * temp_x[k][2] + x_err_comp_sum[k][2];
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

        memcpy(temp_x_err_comp_sum, x_err_comp_sum, objects_count * 3 * sizeof(real));
        memcpy(temp_v_err_comp_sum, v_err_comp_sum, objects_count * 3 * sizeof(real));
        for (int k = 0; k < objects_count; k++)
        {
            temp_v_err_comp_sum[k][0] += *dt * temp_v[k][0];
            temp_v_err_comp_sum[k][1] += *dt * temp_v[k][1];
            temp_v_err_comp_sum[k][2] += *dt * temp_v[k][2];
            temp_x_err_comp_sum[k][0] += *dt * temp_x[k][0];
            temp_x_err_comp_sum[k][1] += *dt * temp_x[k][1];
            temp_x_err_comp_sum[k][2] += *dt * temp_x[k][2];

            v_1[k][0] = v[k][0] + temp_v_err_comp_sum[k][0];
            v_1[k][1] = v[k][1] + temp_v_err_comp_sum[k][1];
            v_1[k][2] = v[k][2] + temp_v_err_comp_sum[k][2];
            x_1[k][0] = x[k][0] + temp_x_err_comp_sum[k][0];
            x_1[k][1] = x[k][1] + temp_x_err_comp_sum[k][1];
            x_1[k][2] = x[k][2] + temp_x_err_comp_sum[k][2];

            temp_v_err_comp_sum[k][0] += v[k][0] - v_1[k][0];
            temp_v_err_comp_sum[k][1] += v[k][1] - v_1[k][1];
            temp_v_err_comp_sum[k][2] += v[k][2] - v_1[k][2];
            temp_x_err_comp_sum[k][0] += x[k][0] - x_1[k][0];
            temp_x_err_comp_sum[k][1] += x[k][1] - x_1[k][1];
            temp_x_err_comp_sum[k][2] += x[k][2] - x_1[k][2];
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

        if (error <= 1 || *dt == tf * 1e-12)
        {
            // Advance step
            *t += *dt; 
            memcpy(x, x_1, objects_count * 3 * sizeof(real));
            memcpy(v, v_1, objects_count * 3 * sizeof(real));
            *count += 1;

            memcpy(x_err_comp_sum, temp_x_err_comp_sum, objects_count * 3 * sizeof(real));
            memcpy(v_err_comp_sum, temp_v_err_comp_sum, objects_count * 3 * sizeof(real));

            // Store step
            if ((*count + 1) % store_every_n == 0)
            {
                sol_time[*store_count + 1] = *t;
                sol_dt[*store_count + 1] = *dt;
                for (int j = 0; j < objects_count; j++)
                {
                    sol_state[*store_count + 1][j * 3] = x[j][0];
                    sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                    sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
                }  
                *store_count += 1;
            }

            if (*t == tf)
            {
                sol_time[*store_count] = *t;
                sol_dt[*store_count] = *dt;
                for (int j = 0; j < objects_count; j++)
                {
                    sol_state[*store_count][j * 3] = x[j][0];
                    sol_state[*store_count][j * 3 + 1] = x[j][1];
                    sol_state[*store_count][j * 3 + 2] = x[j][2];
                    sol_state[*store_count][objects_count * 3 + j * 3] = v[j][0];
                    sol_state[*store_count][objects_count * 3 + j * 3 + 1] = v[j][1];
                    sol_state[*store_count][objects_count * 3 + j * 3 + 2] = v[j][2];
                }  
            }

            // Check buffer size and quit if full
            if ((*store_count + 1) == len_sol_time)
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
                free(temp_x_err_comp_sum);
                free(temp_v_err_comp_sum);
                return 1;
            } 
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

        if (dt_new / tf < 1e-12)
        {
            *dt = tf * 1e-12;
        }

        // Correct overshooting
        if (*t + *dt > tf)
        {
            *dt = tf - *t;
        }
        
        // Exit if finished
        if (*t >= tf)
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
            free(temp_x_err_comp_sum);
            free(temp_v_err_comp_sum);
            return 2;
        }

        // Exit to update progress bar
        if ((int) (*t / tf * 100.0) > progress_percentage)
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
            free(temp_x_err_comp_sum);
            free(temp_v_err_comp_sum);
            return 0;
        }
    }
}

WIN32DLL_API int ias15(
    int objects_count, 
    int dim_nodes,
    const real *restrict nodes,
    const real *restrict aux_c, 
    const real *restrict aux_r,
    real *restrict aux_b0,
    real *restrict aux_b,
    real *restrict aux_g,
    real *restrict aux_e,
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real (*restrict a)[3], 
    const real *restrict m, 
    real G,     
    real *restrict t, 
    real *restrict dt, 
    real tf, 
    int store_every_n,
    int *store_count,
    int *count, 
    real tolerance,
    real tolerance_pc,
    real (*restrict sol_state)[6 * objects_count],
    int len_sol_time,
    real *restrict sol_time,
    real *restrict sol_dt,
    real safety_fac,
    real exponent,
    int *restrict ias15_refine_flag,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3]
)
{
    // Round down current progress percentage as int
    int progress_percentage = (int) (*t / tf * 100.0);

    int dim_nodes_minus_1 = dim_nodes - 1;
    int dim_nodes_minus_2 = dim_nodes - 2;

    // Arrays for ias15_step
    real *aux_a = malloc(dim_nodes * objects_count * 3 * sizeof(real));
    real (*x_step)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*v_step)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a_step)[3] = malloc(objects_count * 3 * sizeof(real));
    real *delta_b7 = malloc(objects_count * 3 * sizeof(real));

    // Array for compute aux_g
    real *F = malloc(8 * objects_count * 3 * sizeof(real));

    // Array for refine aux_b
    real *delta_aux_b = malloc(dim_nodes_minus_1 * objects_count * 3 * sizeof(real));

    // Arrays for compensated summation
    real (*temp_x_err_comp_sum)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_v_err_comp_sum)[3] = malloc(objects_count * 3 * sizeof(real));

    while (1)
    {
        ias15_step(
            objects_count,
            dim_nodes,
            dim_nodes_minus_1,
            dim_nodes_minus_2,
            x,
            v,
            a,
            m,
            G,
            t,
            dt,
            tf,
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
            delta_aux_b,
            x_err_comp_sum, 
            v_err_comp_sum,
            temp_x_err_comp_sum,
            temp_v_err_comp_sum
        );

        // Update count
        *count += 1;

            // Store step
            if ((*count + 1) % store_every_n == 0)
            {
                sol_time[*store_count + 1] = *t;
                sol_dt[*store_count + 1] = *dt;
                for (int j = 0; j < objects_count; j++)
                {
                    sol_state[*store_count + 1][j * 3] = x[j][0];
                    sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                    sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
                }  
                *store_count += 1;
            }

            if (*t == tf)
            {
                sol_time[*store_count] = *t;
                sol_dt[*store_count] = *dt;
                for (int j = 0; j < objects_count; j++)
                {
                    sol_state[*store_count][j * 3] = x[j][0];
                    sol_state[*store_count][j * 3 + 1] = x[j][1];
                    sol_state[*store_count][j * 3 + 2] = x[j][2];
                    sol_state[*store_count][objects_count * 3 + j * 3] = v[j][0];
                    sol_state[*store_count][objects_count * 3 + j * 3 + 1] = v[j][1];
                    sol_state[*store_count][objects_count * 3 + j * 3 + 2] = v[j][2];
                }  
            }

        // Check buffer size and quit if full
        if ((*store_count + 1) == len_sol_time)
        {   
            free(aux_a);
            free(x_step);
            free(v_step);
            free(a_step);
            free(delta_b7); 
            free(F);
            free(delta_aux_b);
            free(temp_x_err_comp_sum);
            free(temp_v_err_comp_sum);
            return 1;
        } 

        // End simulation as t = tf
        if (*t >= tf)
        {
            free(aux_a);
            free(x_step);
            free(v_step);
            free(a_step);
            free(delta_b7); 
            free(F);
            free(delta_aux_b);
            free(temp_x_err_comp_sum);
            free(temp_v_err_comp_sum);
            return 2;
        }

        // Exit to update progress bar
        if ((int) (*t / tf * 100.0) > progress_percentage)
        {
            free(aux_a);
            free(x_step);
            free(v_step);
            free(a_step);
            free(delta_b7); 
            free(F);
            free(delta_aux_b);
            free(temp_x_err_comp_sum);
            free(temp_v_err_comp_sum);
            return 0;
        }
    }
}

// Advance IAS15 for one step
WIN32DLL_API void ias15_step(
    int objects_count,
    int dim_nodes,
    int dim_nodes_minus_1,
    int dim_nodes_minus_2,
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    const real *restrict m,
    real G,
    real *restrict t,
    real *restrict dt,
    real tf,
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
    real *restrict delta_aux_b,
    real (*restrict x_err_comp_sum)[3], 
    real (*restrict v_err_comp_sum)[3],
    real (*restrict temp_x_err_comp_sum)[3], 
    real (*restrict temp_v_err_comp_sum)[3]
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
                ias15_approx_pos_aux(objects_count, x, x0, v0, a0, nodes[i], aux_b, *dt, x_err_comp_sum);
                ias15_approx_vel_aux(objects_count, v, v0, a0, nodes[i], aux_b, *dt, v_err_comp_sum);

                // Evaluate force function and store result
                acceleration(objects_count, x, a, m, G);
                memcpy(&aux_a[i * objects_count * 3], a, objects_count * 3 * sizeof(real));
                
                ias15_compute_aux_g(objects_count, dim_nodes, aux_g, aux_r, aux_a, i, F);
                ias15_compute_aux_b(objects_count, dim_nodes_minus_1, aux_b, aux_g, aux_c, i);
            }

            // Estimate convergence
            for (int i = 0; i < objects_count; i++)
            {
                delta_b7[i * 3 + 0] = aux_b[dim_nodes_minus_2 * objects_count * 3 + i * 3 + 0] - aux_b0[dim_nodes_minus_2 * objects_count * 3 + i * 3 + 0];
                delta_b7[i * 3 + 1] = aux_b[dim_nodes_minus_2 * objects_count * 3 + i * 3 + 1] - aux_b0[dim_nodes_minus_2 * objects_count * 3 + i * 3 + 1];
                delta_b7[i * 3 + 2] = aux_b[dim_nodes_minus_2 * objects_count * 3 + i * 3 + 2] - aux_b0[dim_nodes_minus_2 * objects_count * 3 + i * 3 + 2];
            }
            memcpy(aux_b0, aux_b, dim_nodes_minus_1 * objects_count * 3 * sizeof(real));
            if ((abs_max_vec(delta_b7, objects_count * 3) / abs_max_vec(&aux_a[dim_nodes_minus_1 * objects_count * 3], objects_count * 3)) < tolerance_pc)
            {
                break;
            }
        }

        // Advance step
        memcpy(temp_x_err_comp_sum, x_err_comp_sum, objects_count * 3 * sizeof(real));
        memcpy(temp_v_err_comp_sum, v_err_comp_sum, objects_count * 3 * sizeof(real));

        ias15_approx_pos_step(objects_count, x, x0, v0, a0, aux_b, *dt, temp_x_err_comp_sum);
        ias15_approx_vel_step(objects_count, v, v0, a0, aux_b, *dt, temp_v_err_comp_sum);
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
        if (error <= 1 || *dt == tf * 1e-12)
        {
            // Report accepted step
            ias15_integrate_flag = 1;
            *t += *dt;

            ias15_refine_aux_b(objects_count, dim_nodes_minus_1, aux_b, aux_e, delta_aux_b, *dt, dt_new, *ias15_refine_flag);
            *ias15_refine_flag = 1;

            memcpy(x_err_comp_sum, temp_x_err_comp_sum, objects_count * 3 * sizeof(real));
            memcpy(v_err_comp_sum, temp_v_err_comp_sum, objects_count * 3 * sizeof(real));

            if (*t >= tf)
            {
                memcpy(x0, x, objects_count * 3 * sizeof(real));
                memcpy(v0, v, objects_count * 3 * sizeof(real));
                memcpy(a0, a, objects_count * 3 * sizeof(real));     
                break;  
            }
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

        if (dt_new / tf < 1e-12)
        {
            *dt = tf * 1e-12;
        }

        // Correct overshooting
        if (*t + *dt > tf)
        {
            *dt = tf - *t;
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

// Calculate the auxiliary position used to calculate aux_b and aux_g
WIN32DLL_API void ias15_approx_pos_aux(
    int objects_count,
    real (*restrict x)[3],
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real node,
    real *restrict aux_b,
    real dt,
    real (*restrict x_err_comp_sum)[3]
)
{   
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            x[j][k] = x0[j][k] + x_err_comp_sum[j][k] + dt * node * (
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

// Calculate the auxiliary velocity used to calculate aux_b and aux_g
WIN32DLL_API void ias15_approx_vel_aux(
    int objects_count,
    real (*restrict v)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real node,
    real *restrict aux_b,
    real dt,
    real (*restrict v_err_comp_sum)[3]
)
{
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            v[j][k] = v0[j][k] + v_err_comp_sum[j][k] + dt * node * (
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

// Calculate the position of the next step
WIN32DLL_API void ias15_approx_pos_step(
    int objects_count,
    real (*restrict x)[3],
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real *restrict aux_b,
    real dt,
    real (*restrict temp_x_err_comp_sum)[3]
)
{   
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            temp_x_err_comp_sum[j][k] += dt * (
                v0[j][k] + dt * (a0[j][k]
                    + aux_b[0 * objects_count * 3 + j * 3 + k] / 3.0
                    + aux_b[1 * objects_count * 3 + j * 3 + k] / 6.0
                    + aux_b[2 * objects_count * 3 + j * 3 + k] / 10.0
                    + aux_b[3 * objects_count * 3 + j * 3 + k] / 15.0
                    + aux_b[4 * objects_count * 3 + j * 3 + k] / 21.0
                    + aux_b[5 * objects_count * 3 + j * 3 + k] / 28.0 
                    + aux_b[6 * objects_count * 3 + j * 3 + k] / 36.0
                )
                / 2.0
            );

            x[j][k] = x0[j][k] + temp_x_err_comp_sum[j][k];
            temp_x_err_comp_sum[j][k] += (x0[j][k] - x[j][k]);
        }
    }
}

// Calculate the velocity of the next step
WIN32DLL_API void ias15_approx_vel_step(
    int objects_count,
    real (*restrict v)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    real *restrict aux_b,
    real dt,
    real (*restrict temp_v_err_comp_sum)[3]
)
{
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            temp_v_err_comp_sum[j][k] += dt * (
                a0[j][k]
                + aux_b[0 * objects_count * 3 + j * 3 + k] / 2.0
                + aux_b[1 * objects_count * 3 + j * 3 + k] / 3.0
                + aux_b[2 * objects_count * 3 + j * 3 + k] / 4.0
                + aux_b[3 * objects_count * 3 + j * 3 + k] / 5.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] / 6.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] / 7.0 
                + aux_b[6 * objects_count * 3 + j * 3 + k] / 8.0
            );
            v[j][k] = v0[j][k] + temp_v_err_comp_sum[j][k];
            temp_v_err_comp_sum[j][k] += (v0[j][k] - v[j][k]);
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
