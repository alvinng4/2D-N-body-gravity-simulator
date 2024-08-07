#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "acceleration.h"

/**
 * \brief IAS15 integrator
 * 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors of all objects
 * \param v Array of velocity vectors of all objects
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param t Pointer to the current simulation time
 * \param tf Total time to be integrated
 * \param input_tolerance Tolerance of the integrator
 * \param acceleration_method Method to calculate acceleration
 * \param barnes_hut_theta Theta parameter for Barnes-Hut algorithm
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
int ias15(
    int objects_count,
    real *restrict x,
    real *restrict v,
    real *restrict m,
    real G,
    double *restrict t,
    double tf, 
    double input_tolerance,
    const char *restrict acceleration_method,
    real barnes_hut_theta,
    int store_every_n,
    int *restrict store_count,
    const int storing_method,
    const char *restrict flush_path,
    Solutions *restrict solution,
    bool *restrict is_exit
);

/**
 * \brief Initialize the radau spacing nodes for IAS15
 * 
 * \param nodes 1D array of size 8 to be modified
 * 
 * \return None
 */
void ias15_radau_spacing(real *restrict nodes);

/**
 * \brief Initialize the auxiliary coefficients aux_c for IAS15
 * 
 * \param aux_c 1D array of length 49 to be modified
 * 
 * \return None
 */
void ias15_aux_c(real *restrict aux_c);

/**
 * \brief Initialize auxiliary coefficients aux_r for IAS15
 * 
 * \param aux_r 1D array of size 64 to be modified
 * 
 * \return None
 */
void ias15_aux_r(real *aux_r);

/**
 * \brief Calculate the initial time step for IAS15 integrator
 * 
 * \param objects_count Number of objects in the system
 * \param power Power of the integrator
 * \param x Array of position vectors of all objects
 * \param v Array of velocity vectors of all objects
 * \param a Array of acceleration vectors of all objects
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param acceleration_method_flag Method to calculate acceleration (int flag)
 * \param barnes_hut_theta Theta parameter for Barnes-Hut algorithm
 * 
 * \return initial dt for IAS15 integrator
 * \retval -1.0 If failed to allocate memory for calculation
 */
real ias15_initial_dt(
    int objects_count,
    int power,
    real *restrict x,
    real *restrict v,
    real *restrict a,
    const real *m,
    real G,
    int acceleration_method_flag,
    real barnes_hut_theta
);

/**
 * \brief Advance IAS15 for one step.
 * 
 * Most of the function parameters are needed in subsequent functions
 * called by ias15_step. They are passed to avoid repeating memory 
 * allocation in order to improve performance.
 * 
 * \param objects_count Number of objects in the system
 * \param dim_nodes Dimension of nodes for the predictor-corrector algorithm
 * \param dim_nodes_minus_1 Dimension of nodes for the predictor-corrector algorithm minus one
 * \param dim_nodes_minus_2 Dimension of nodes for the predictor-corrector algorithm minus two
 * \param x0 Array of position vectors from the last time step
 * \param v0 Array of velocity vectors from the last time step
 * \param a0 Array of acceleration vectors from the last time step
 * \param m Array of masses for all objects
 * \param G Gravitational constant
 * \param t Pointer to the current simulation time
 * \param dt Input of current time step of the system
 * \param dt_old Actual current time step of the system after the function returns
 * \param tf Total time to be integrated
 * \param nodes Radau spacing dodes for the predictor-corrector algorithm
 * \param aux_b0 Array of auxiliary coefficients b0
 * \param aux_b Array of auxiliary coefficients b
 * \param aux_c Array of auxiliary coefficients c
 * \param aux_e Array of auxiliary coefficients e
 * \param aux_g Array of auxiliary coefficients g
 * \param aux_r Array of auxiliary coefficients r 
 * \param tolerance Tolerance of the integrator
 * \param tolerance_pc Tolerance of the predictor-corrector algorithm
 * \param exponent  Exponent of step-size control
 * \param safety_fac Safety factor for step-size control
 * \param ias15_refine_flag Flag for refining the auxiliary coefficients b
 * \param aux_a Array of auxiliary accelerations a
 * \param x Array of position vectors to be modified
 * \param v Array of velocity vectors to be modified
 * \param a Array of acceleration vectors to be modified
 * \param delta_b7 Array of delta b7
 * \param F Helper array for compute_aux_g
 * \param delta_aux_b Helper array for refine_aux_b
 * \param x_err_comp_sum Array of round off errors of position vectors
 * \param v_err_comp_sum Array of round off errors of velocity vectors
 * \param temp_x_err_comp_sum Temporary array of round off errors of position vectors
 * \param temp_v_err_comp_sum Temporary array of round off errors of velocity vectors
 * \param acceleration_method_flag Method to calculate acceleration (int flag)
 * \param barnes_hut_theta Theta parameter for Barnes-Hut algorithm
 * 
 * \return None
 */
void ias15_step(
    int objects_count,
    int dim_nodes,
    int dim_nodes_minus_1,
    int dim_nodes_minus_2,
    real *restrict x0,
    real *restrict v0,
    real *restrict a0,
    const real *restrict m,
    real G,
    real *restrict t,
    real *restrict dt,
    real *restrict dt_old,
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
    real *restrict x,
    real *restrict v,
    real *restrict a,
    real *restrict delta_b7,
    real *restrict F,
    real *restrict delta_aux_b,
    real *restrict x_err_comp_sum, 
    real *restrict v_err_comp_sum,
    real *restrict temp_x_err_comp_sum,
    real *restrict temp_v_err_comp_sum,
    int acceleration_method_flag,
    real barnes_hut_theta
);

/**
 * \brief Calculate position for the predictor-corrector algorithm
 * 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors to be modified
 * \param x0 Array of position vectors from the last time step
 * \param v0 Array of velocity vectors from the last time step
 * \param a0 Array of acceleration vectors from the last time step
 * \param node Current node of the predictor-corrector algorithm
 * \param aux_b Auxiliary b array
 * \param dt Current time step of the system
 * \param x_err_comp_sum Array of round off errors of position vectors 
 *                       for compensated summation
 * 
 * \return None
 */
void ias15_approx_pos_pc(
    int objects_count,
    real *restrict x,
    real *restrict x0,
    real *restrict v0,
    real *restrict a0,
    real node,
    real *restrict aux_b,
    real dt,
    real *restrict x_err_comp_sum
);

/**
 * \brief Calculate velocity for the predictor-corrector algorithm
 *
 * \param objects_count Number of objects in the system
 * \param v Array of velocity vectors to be modified
 * \param v0 Array of velocity vectors from the last time step
 * \param a0 Array of acceleration vectors from the last time step
 * \param node Current node of the predictor-corrector algorithm
 * \param aux_b Auxiliary b array
 * \param dt Current time step of the system
 * \param v_err_comp_sum Array of round off errors of velocity vectors 
 *                       for compensated summation
 * 
 * \return None
 */
void ias15_approx_vel_pc(
    int objects_count,
    real *restrict v,
    real *restrict v0,
    real *restrict a0,
    real node,
    real *restrict aux_b,
    real dt,
    real *restrict v_err_comp_sum
);

/**
 * \brief Calculate position for the next time step
 * 
 * \param objects_count Number of objects in the system
 * \param x Array of position vectors to be modified
 * \param x0 Array of position vectors from the last time step
 * \param v0 Array of velocity vectors from the last time step
 * \param a0 Array of acceleration vectors from the last time step
 * \param aux_b Auxiliary coefficients b
 * \param dt Current time step of the system
 * \param temp_x_err_comp_sum Temporary array of round off errors of position vectors 
 *                            for compensated summation. This array will be
 *                            modified.
 * 
 * \return None
 */
void ias15_approx_pos_step(
    int objects_count,
    real *restrict x,
    real *restrict x0,
    real *restrict v0,
    real *restrict a0,
    real *restrict aux_b,
    real dt,
    real *restrict temp_x_err_comp_sum
);

/**
 * \brief Calculate velocity for the next time step
 * 
 * \param objects_count Number of objects in the system
 * \param v Array of velocity vectors to be modified
 * \param v0 Array of velocity vectors from the last time step
 * \param a0 Array of acceleration vectors from the last time step
 * \param aux_b Auxiliary coefficients b
 * \param dt Current time step of the system
 * \param temp_v_err_comp_sum Temporary array of round off errors of velocity vectors 
 *                            for compensated summation. This array will be
 *                            modified.
 * 
 * \return None
 */
void ias15_approx_vel_step(
    int objects_count,
    real *restrict v,
    real *restrict v0,
    real *restrict a0,
    real *restrict aux_b,
    real dt,
    real *restrict temp_v_err_comp_sum
);

/**
 * \brief Calculate the auxiliary coefficients b for IAS15 integrator
 * 
 * \param objects_count Number of objects in the system
 * \param dim_nodes_minus_1 Dimension of nodes minus one
 * \param aux_b Array of auxiliary coefficients b to be modified
 * \param aux_g Array of auxiliary coefficients g
 * \param aux_c Array of auxiliary coefficients c
 * \param i Current iteration of nodes of the predictor-corrector algorithm
 */
void ias15_compute_aux_b(
    int objects_count,
    int dim_nodes_minus_1,
    real *restrict aux_b,
    const real *restrict aux_g,
    const real *restrict aux_c,
    int i
);

/**
 * \brief Calculate the auxiliary coefficients g for IAS15 integrator
 * 
 * \param objects_count Number of objects in the system
 * \param dim_nodes Dimension of nodes
 * \param aux_g Array of auxiliary coefficients g to be modified
 * \param aux_r Array of auxiliary coefficients r
 * \param aux_a Array of auxiliary accelerations a
 * \param i Current iteration of nodes of the predictor-corrector algorithm
 * \param F Helper array for this function
 */
void ias15_compute_aux_g(
    int objects_count,
    int dim_nodes,
    real *restrict aux_g,
    const real *restrict aux_r,
    const real *restrict aux_a,
    int i,
    real *restrict F
);

/**
 * \brief Refine the auxiliary coefficients b for IAS15 integrator
 * 
 * \param objects_count Number of objects in the system
 * \param dim_nodes_minus_1 Dimension of nodes minus one
 * \param aux_b Array of auxiliary coefficients b to be modified
 * \param aux_e Array of auxiliary coefficients e
 * \param delta_aux_b Helper array for this function
 * \param dt Current time step of the system
 * \param dt_new Next Time Step of the system
 * \param ias15_refine_flag Helper flag for this function
 * 
 * \return None
 */
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


WIN32DLL_API int ias15(
    int objects_count,
    real *restrict x,
    real *restrict v,
    real *restrict m,
    real G,
    double *restrict t,
    double tf, 
    double input_tolerance,
    const char *restrict acceleration_method,
    real barnes_hut_theta,
    int store_every_n,
    int *restrict store_count,
    const int storing_method,
    const char *restrict flush_path,
    Solutions *restrict solution,
    bool *restrict is_exit
)
{   
    int acceleration_method_flag;
    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        acceleration_method_flag = 0;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        acceleration_method_flag = 1;
    }
    else if (strcmp(acceleration_method, "barnes-hut") == 0)
    {
        acceleration_method_flag = 2;
    }
    else
    {
        fprintf(stderr, "Error: acceleration method not recognized\n");
        goto err_acc_method;
    }
    
    real *a = malloc(objects_count * 3 * sizeof(real));
    if (!a)
    {
        fprintf(stderr, "Error: Failed to allocate memory for acceleration vectors\n");
        goto err_a;
    }

    // Recommended tolerance: 1e-9
    real tolerance = input_tolerance;

    // Safety factors for step-size control
    real safety_fac = 0.25;

    // For fixed step integration, choose exponent = 0
    real exponent = 1.0 / 7.0;

    // Tolerance of predictor-corrector algorithm
    real tolerance_pc = 1e-16;
    
    // Initialize auxiliary variables
    int dim_nodes = 8;
    int dim_nodes_minus_1 = 7;
    int dim_nodes_minus_2 = 6;
    real *restrict nodes = calloc(dim_nodes, sizeof(real));
    real *restrict aux_c = calloc(7 * 7, sizeof(real));
    real *restrict aux_r = calloc(8 * 8, sizeof(real));
    real *restrict aux_b0 = calloc((dim_nodes - 1) * objects_count * 3, sizeof(real));
    real *restrict aux_b = calloc((dim_nodes - 1) * objects_count * 3, sizeof(real));
    real *restrict aux_g = calloc((dim_nodes - 1) * objects_count * 3, sizeof(real));
    real *restrict aux_e = calloc((dim_nodes - 1) * objects_count * 3, sizeof(real));

    if (!nodes || !aux_c || !aux_r || !aux_b0 || !aux_b || !aux_g || !aux_e)
    {
        fprintf(stderr, "Error: Failed to allocate memory for auxiliary variables\n");
        goto err_aux_memory;
    }
    ias15_radau_spacing(nodes);
    ias15_aux_c(aux_c);
    ias15_aux_r(aux_r);

    int ias15_refine_flag = 0;

    // Arrays for ias15_step
    real *restrict aux_a = calloc(dim_nodes * objects_count * 3, sizeof(real));
    real *restrict x_step = calloc(objects_count * 3, sizeof(real));
    real *restrict v_step = calloc(objects_count * 3, sizeof(real));
    real *restrict a_step = calloc(objects_count * 3, sizeof(real));
    real *restrict delta_b7 = calloc(objects_count * 3, sizeof(real));

    // Array for compute aux_g
    real *restrict F = calloc(8 * objects_count * 3, sizeof(real));

    // Array for refine aux_b
    real *restrict delta_aux_b = calloc(dim_nodes_minus_1 * objects_count * 3, sizeof(real));

    // Arrays for compensated summation
    real *restrict x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict v_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict temp_x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *restrict temp_v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (
        !aux_a || 
        !x_step ||
        !v_step ||
        !a_step ||
        !delta_b7 ||
        !F ||
        !delta_aux_b ||
        !x_err_comp_sum ||
        !v_err_comp_sum ||
        !temp_x_err_comp_sum ||
        !temp_v_err_comp_sum
    )
    {
        fprintf(stderr, "Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    int64 count = 1;    // 1 for t0
    real dt_old = ias15_initial_dt(objects_count, 15, x, v, a, m, G, acceleration_method_flag, barnes_hut_theta);
    real dt = dt_old;

    FILE *flush_file = NULL;
    double *sol_state = NULL;
    double *sol_time = NULL;
    double *sol_dt = NULL;
    int buffer_size = BUFFER_SIZE;

    // For realloc solution output
    double *temp_sol_state = NULL;
    double *temp_sol_time = NULL;
    double *temp_sol_dt = NULL;
    
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
        sol_state = malloc(buffer_size * objects_count * 6 * sizeof(double));
        sol_time = malloc(buffer_size * sizeof(double));
        sol_dt = malloc(buffer_size * sizeof(double));

        if (!sol_state || !sol_time || !sol_dt)
        {
            fprintf(stderr, "Error: Failed to allocate memory for solution output\n");
            goto err_sol_output_memory;
        }

        // Initial value
        memcpy(&sol_state[0], x, objects_count * 3 * sizeof(double));
        memcpy(&sol_state[objects_count * 3], v, objects_count * 3 * sizeof(double));
        sol_time[0] = 0.0;
        if (dt == -1.0)
        {
            goto err_initial_dt_memory;
        }
        sol_dt[0] = dt;
    }

    while (*t < tf)
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
            &dt,
            &dt_old,
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
            &ias15_refine_flag,
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
            temp_v_err_comp_sum,
            acceleration_method_flag,
            barnes_hut_theta
        );

        // Store step
        if (count % store_every_n == 0)
        {
            if (storing_method == 1)
            {
                write_to_csv_file(flush_file, *t, dt, objects_count, x, v, m, G);
            }
            else if (storing_method == 0)
            {
                sol_time[*store_count] = *t;
                sol_dt[*store_count] = dt_old;
                memcpy(&sol_state[*store_count * objects_count * 6], x, objects_count * 3 * sizeof(double));
                memcpy(&sol_state[*store_count * objects_count * 6 + objects_count * 3], v, objects_count * 3 * sizeof(double));
            }
            (*store_count)++;

            // Check buffer size and extend if full
            if ((storing_method == 0) && (*store_count == buffer_size) && (*t < tf))
            {   
                buffer_size *= 2;
                temp_sol_state = realloc(sol_state, buffer_size * objects_count * 6 * sizeof(real));
                temp_sol_time = realloc(sol_time, buffer_size * sizeof(real));
                temp_sol_dt = realloc(sol_dt, buffer_size * sizeof(real));

                if (!temp_sol_state || !temp_sol_time || !temp_sol_dt)
                {
                    fprintf(stderr, "Error: Failed to allocate extra memory to extend array for solution output\n");
                    goto err_sol_output_memory;
                }
                
                sol_state = temp_sol_state;
                sol_time = temp_sol_time;
                sol_dt = temp_sol_dt;
            }
        }
        count++;
    
        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            goto err_user_exit;
        }   
    }

    free(a);
    free(nodes);
    free(aux_c);
    free(aux_r);
    free(aux_b0);
    free(aux_b);
    free(aux_g);
    free(aux_e);
    free(aux_a);
    free(x_step);
    free(v_step);
    free(a_step);
    free(delta_b7); 
    free(F);
    free(delta_aux_b);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    free(temp_x_err_comp_sum);
    free(temp_v_err_comp_sum);
    
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
err_initial_dt_memory:
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
    free(aux_a);
    free(x_step);
    free(v_step);
    free(a_step);
    free(delta_b7);
    free(F);
    free(delta_aux_b);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    free(temp_x_err_comp_sum);
    free(temp_v_err_comp_sum);
err_aux_memory:
    free(nodes);
    free(aux_c);
    free(aux_r);
    free(aux_b0);
    free(aux_b);
    free(aux_g);
    free(aux_e);
err_a:
    free(a);
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

WIN32DLL_API void ias15_radau_spacing(real *restrict nodes)
{
    nodes[0] = 0.0L;
    nodes[1] = 0.056262560536922146465652191032L;
    nodes[2] = 0.180240691736892364987579942809L;
    nodes[3] = 0.352624717113169637373907770171L;
    nodes[4] = 0.547153626330555383001448557652L;
    nodes[5] = 0.734210177215410531523210608306L;
    nodes[6] = 0.885320946839095768090359762932L;
    nodes[7] = 0.977520613561287501891174500429L;
}

WIN32DLL_API void ias15_aux_c(real *restrict aux_c)
{
    for (int i = 0; i < 7; i++)
    {
        aux_c[i * 7 + i] = 1.0L;
    }

    aux_c[1 * 7 + 0] = -0.0562625605369221464656522L;

    aux_c[2 * 7 + 0] = 0.01014080283006362998648180399549641417413495311078L;
    aux_c[2 * 7 + 1] = -0.2365032522738145114532321L;

    aux_c[3 * 7 + 0] = -0.0035758977292516175949344589284567187362040464593728L;
    aux_c[3 * 7 + 1] = 0.09353769525946206589574845561035371499343547051116L;
    aux_c[3 * 7 + 2] = -0.5891279693869841488271399L;

    aux_c[4 * 7 + 0] = 0.0019565654099472210769005672379668610648179838140913L;
    aux_c[4 * 7 + 1] = -0.054755386889068686440808430671055022602028382584495L;
    aux_c[4 * 7 + 2] = 0.41588120008230686168862193041156933067050816537030L;
    aux_c[4 * 7 + 3] = -1.1362815957175395318285885L;

    aux_c[5 * 7 + 0] = -0.0014365302363708915424459554194153247134438571962198L;
    aux_c[5 * 7 + 1] = 0.042158527721268707707297347813203202980228135395858L;
    aux_c[5 * 7 + 2] = -0.36009959650205681228976647408968845289781580280782L;
    aux_c[5 * 7 + 3] = 1.2501507118406910258505441186857527694077565516084L;
    aux_c[5 * 7 + 4] = -1.8704917729329500633517991L;

    aux_c[6 * 7 + 0] = 0.0012717903090268677492943117622964220889484666147501L;
    aux_c[6 * 7 + 1] = -0.038760357915906770369904626849901899108502158354383L;
    aux_c[6 * 7 + 2] = 0.36096224345284598322533983078129066420907893718190L;
    aux_c[6 * 7 + 3] = -1.4668842084004269643701553461378480148761655599754L;
    aux_c[6 * 7 + 4] = 2.9061362593084293014237914371173946705384212479246L;
    aux_c[6 * 7 + 5] = -2.7558127197720458314421589L;
}

WIN32DLL_API void ias15_aux_r(real *aux_r)
{
    aux_r[1 * 8 + 0] = 17.773808914078000840752659565672904106978971632681L;
    aux_r[2 * 8 + 0] = 5.5481367185372165056928216140765061758579336941398L;
    aux_r[3 * 8 + 0] = 2.8358760786444386782520104428042437400879003147949L;
    aux_r[4 * 8 + 0] = 1.8276402675175978297946077587371204385651628457154L;
    aux_r[5 * 8 + 0] = 1.3620078160624694969370006292445650994197371928318L;
    aux_r[6 * 8 + 0] = 1.1295338753367899027322861542728593509768148769105L;
    aux_r[7 * 8 + 0] = 1.0229963298234867458386119071939636779024159134103L;

    aux_r[2 * 8 + 1] = 8.0659386483818866885371256689687154412267416180207L;
    aux_r[3 * 8 + 1] = 3.3742499769626352599420358188267460448330087696743L;
    aux_r[4 * 8 + 1] = 2.0371118353585847827949159161566554921841792590404L;
    aux_r[5 * 8 + 1] = 1.4750402175604115479218482480167404024740127431358L;
    aux_r[6 * 8 + 1] = 1.2061876660584456166252036299646227791474203527801L;
    aux_r[7 * 8 + 1] = 1.0854721939386423840467243172568913862030118679827L;

    aux_r[3 * 8 + 2] = 5.8010015592640614823286778893918880155743979164251L;
    aux_r[4 * 8 + 2] = 2.7254422118082262837742722003491334729711450288807L;
    aux_r[5 * 8 + 2] = 1.8051535801402512604391147435448679586574414080693L;
    aux_r[6 * 8 + 2] = 1.4182782637347391537713783674858328433713640692518L;
    aux_r[7 * 8 + 2] = 1.2542646222818777659905422465868249586862369725826L;

    aux_r[4 * 8 + 3] = 5.1406241058109342286363199091504437929335189668304L;
    aux_r[5 * 8 + 3] = 2.6206449263870350811541816031933074696730227729812L;
    aux_r[6 * 8 + 3] = 1.8772424961868100972169920283109658335427446084411L;
    aux_r[7 * 8 + 3] = 1.6002665494908162609916716949161150366323259154408L;

    aux_r[5 * 8 + 4] = 5.3459768998711075141214909632277898045770336660354L;
    aux_r[6 * 8 + 4] = 2.9571160172904557478071040204245556508352776929762L;
    aux_r[7 * 8 + 4] = 2.3235983002196942228325345451091668073608955835034L;

    aux_r[6 * 8 + 5] = 6.6176620137024244874471284891193925737033291491748L;
    aux_r[7 * 8 + 5] = 4.1099757783445590862385761824068782144723082633980L;

    aux_r[7 * 8 + 6] = 10.846026190236844684706431007823415424143683137181L;
}

WIN32DLL_API real ias15_initial_dt(
    int objects_count,
    int power,
    real *restrict x,
    real *restrict v,
    real *restrict a,
    const real *m,
    real G,
    int acceleration_method_flag,
    real barnes_hut_theta
)
{
    acceleration(acceleration_method_flag, objects_count, x, v, a, m, G, barnes_hut_theta);

    real d_0 = abs_max_vec(x, objects_count * 3);
    real d_1 = abs_max_vec(a, objects_count * 3);
    real d_2;
    real dt_0;
    real dt_1;
    real *x_1 = malloc(objects_count * 3 * sizeof(real));
    real *a_1 = malloc(objects_count * 3 * sizeof(real));

    if (!x_1 || !a_1)
    {
        fprintf(stderr, "Error: Failed to allocate memory for ias15_initial_dt()\n");
        free(x_1);
        free(a_1);
        return -1.0;
    }

    if (d_0 < 1e-5 || d_1 < 1e-5)
    {
        dt_0 = 1e-6;
    }
    else
    {
        dt_0 = 0.01 * (d_0 / d_1);
    }
    
    for (int i = 0; i < objects_count; i++)
    {
        x_1[i * 3 + 0] = x[i * 3 + 0] + dt_0 * v[i * 3 + 0];
        x_1[i * 3 + 1] = x[i * 3 + 1] + dt_0 * v[i * 3 + 1];
        x_1[i * 3 + 2] = x[i * 3 + 2] + dt_0 * v[i * 3 + 2];
    }
    acceleration(acceleration_method_flag, objects_count, x_1, v, a_1, m, G, barnes_hut_theta);

    for (int i = 0; i < objects_count; i++)
    {
        a_1[i * 3 + 0] = a_1[i * 3 + 0] - a[i * 3 + 0];
        a_1[i * 3 + 1] = a_1[i * 3 + 1] - a[i * 3 + 1];
        a_1[i * 3 + 2] = a_1[i * 3 + 2] - a[i * 3 + 2];
    }
    d_2 = abs_max_vec(a_1, objects_count * 3) / dt_0;

    if (fmax(d_1, d_2) <= 1e-15)
    {
        dt_1 = fmax(1e-6, dt_0 * 1e-3);
    }
    else
    {
        dt_1 = pow((0.01 / fmax(d_1, d_2)), 1.0 / (1.0 + power));
    }
    
    free(x_1);
    free(a_1);

    return fmin(100.0 * dt_0, dt_1);
}

WIN32DLL_API void ias15_step(
    int objects_count,
    int dim_nodes,
    int dim_nodes_minus_1,
    int dim_nodes_minus_2,
    real *restrict x0,
    real *restrict v0,
    real *restrict a0,
    const real *restrict m,
    real G,
    real *restrict t,
    real *restrict dt,
    real *restrict dt_old,
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
    real *restrict x,
    real *restrict v,
    real *restrict a,
    real *restrict delta_b7,
    real *restrict F,
    real *restrict delta_aux_b,
    real *restrict x_err_comp_sum, 
    real *restrict v_err_comp_sum,
    real *restrict temp_x_err_comp_sum,
    real *restrict temp_v_err_comp_sum,
    int acceleration_method_flag,
    real barnes_hut_theta
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
                ias15_approx_pos_pc(objects_count, x, x0, v0, a0, nodes[i], aux_b, *dt, x_err_comp_sum);
                ias15_approx_vel_pc(objects_count, v, v0, a0, nodes[i], aux_b, *dt, v_err_comp_sum);

                // Evaluate force function and store result
                acceleration(acceleration_method_flag, objects_count, x, v, &aux_a[i * objects_count * 3], m, G, barnes_hut_theta);

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
        memcpy(temp_x_err_comp_sum, x_err_comp_sum, objects_count * 3 * sizeof(real));
        memcpy(temp_v_err_comp_sum, v_err_comp_sum, objects_count * 3 * sizeof(real));

        ias15_approx_pos_step(objects_count, x, x0, v0, a0, aux_b, *dt, temp_x_err_comp_sum);
        ias15_approx_vel_step(objects_count, v, v0, a0, aux_b, *dt, temp_v_err_comp_sum);
        acceleration(acceleration_method_flag, objects_count, x, v, a, m, G, barnes_hut_theta);

        // Estimate relative error
        error_b7 = abs_max_vec(&aux_b[dim_nodes_minus_2 * objects_count * 3], objects_count * 3) / abs_max_vec(a, objects_count * 3);
        error = pow((error_b7 / tolerance), exponent);

        // Step size of the next step for refine aux b
        *dt_old = *dt;
        if (error < 1e-10)
        {
            error = 1e-10; // Prevent error from being too small
        }
        dt_new = *dt / error;

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
        }

        // Actual step size for the next step
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

WIN32DLL_API void ias15_approx_pos_pc(
    int objects_count,
    real *restrict x,
    real *restrict x0,
    real *restrict v0,
    real *restrict a0,
    real node,
    real *restrict aux_b,
    real dt,
    real *restrict x_err_comp_sum
)
{   
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            /*
            *   Warning: Combining both statements would increase floating point error
            *            e.g. x[j][k] = x0[j][k] + x_err_comp_sum[j][k] + ...   (WRONG)
            */

            x[j * 3 + k] = x0[j * 3 + k];     
            x[j * 3 + k] += x_err_comp_sum[j * 3 + k] + dt * node * (
                v0[j * 3 + k]
                + dt
                * node
                * (
                    a0[j * 3 + k]
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

WIN32DLL_API void ias15_approx_vel_pc(
    int objects_count,
    real *restrict v,
    real *restrict v0,
    real *restrict a0,
    real node,
    real *restrict aux_b,
    real dt,
    real *restrict v_err_comp_sum
)
{
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            /*
            *   Warning: Combining both statements would increase floating point error
            *            e.g. v[j][k] = v0[j][k] + v_err_comp_sum[j][k] + ...   (WRONG)
            */

            v[j * 3 + k] = v0[j * 3 + k];
            v[j * 3 + k] += v_err_comp_sum[j * 3 + k] + dt * node * (
                a0[j * 3 + k]
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

WIN32DLL_API void ias15_approx_pos_step(
    int objects_count,
    real *restrict x,
    real *restrict x0,
    real *restrict v0,
    real *restrict a0,
    real *restrict aux_b,
    real dt,
    real *restrict temp_x_err_comp_sum
)
{   
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            temp_x_err_comp_sum[j * 3 + k] += dt * (
                v0[j * 3 + k] + dt * (a0[j * 3 + k]
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

            x[j * 3 + k] = x0[j * 3 + k] + temp_x_err_comp_sum[j * 3 + k];
            temp_x_err_comp_sum[j * 3 + k] += (x0[j * 3 + k] - x[j * 3 + k]);
        }
    }
}

WIN32DLL_API void ias15_approx_vel_step(
    int objects_count,
    real *restrict v,
    real *restrict v0,
    real *restrict a0,
    real *restrict aux_b,
    real dt,
    real *restrict temp_v_err_comp_sum
)
{
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            temp_v_err_comp_sum[j * 3 + k] += dt * (
                a0[j * 3 + k]
                + aux_b[0 * objects_count * 3 + j * 3 + k] / 2.0
                + aux_b[1 * objects_count * 3 + j * 3 + k] / 3.0
                + aux_b[2 * objects_count * 3 + j * 3 + k] / 4.0
                + aux_b[3 * objects_count * 3 + j * 3 + k] / 5.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] / 6.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] / 7.0 
                + aux_b[6 * objects_count * 3 + j * 3 + k] / 8.0
            );
            v[j * 3 + k] = v0[j * 3 + k] + temp_v_err_comp_sum[j * 3 + k];
            temp_v_err_comp_sum[j * 3 + k] += (v0[j * 3 + k] - v[j * 3 + k]);
        }
    }
}

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
        if (i >= 1) 
        {
            aux_b[0 * objects_count * 3 + j * 3 + 0] = (
                aux_c[0 * dim_nodes_minus_1 + 0] * aux_g[0 * objects_count * 3 + j * 3 + 0]
                + aux_c[1 * dim_nodes_minus_1 + 0] * aux_g[1 * objects_count * 3 + j * 3 + 0]
                + aux_c[2 * dim_nodes_minus_1 + 0] * aux_g[2 * objects_count * 3 + j * 3 + 0]
                + aux_c[3 * dim_nodes_minus_1 + 0] * aux_g[3 * objects_count * 3 + j * 3 + 0]
                + aux_c[4 * dim_nodes_minus_1 + 0] * aux_g[4 * objects_count * 3 + j * 3 + 0]
                + aux_c[5 * dim_nodes_minus_1 + 0] * aux_g[5 * objects_count * 3 + j * 3 + 0]
                + aux_c[6 * dim_nodes_minus_1 + 0] * aux_g[6 * objects_count * 3 + j * 3 + 0]
            );
            aux_b[0 * objects_count * 3 + j * 3 + 1] = (
                aux_c[0 * dim_nodes_minus_1 + 0] * aux_g[0 * objects_count * 3 + j * 3 + 1]
                + aux_c[1 * dim_nodes_minus_1 + 0] * aux_g[1 * objects_count * 3 + j * 3 + 1]
                + aux_c[2 * dim_nodes_minus_1 + 0] * aux_g[2 * objects_count * 3 + j * 3 + 1]
                + aux_c[3 * dim_nodes_minus_1 + 0] * aux_g[3 * objects_count * 3 + j * 3 + 1]
                + aux_c[4 * dim_nodes_minus_1 + 0] * aux_g[4 * objects_count * 3 + j * 3 + 1]
                + aux_c[5 * dim_nodes_minus_1 + 0] * aux_g[5 * objects_count * 3 + j * 3 + 1]
                + aux_c[6 * dim_nodes_minus_1 + 0] * aux_g[6 * objects_count * 3 + j * 3 + 1]
            );
            aux_b[0 * objects_count * 3 + j * 3 + 2] = (
                aux_c[0 * dim_nodes_minus_1 + 0] * aux_g[0 * objects_count * 3 + j * 3 + 2]
                + aux_c[1 * dim_nodes_minus_1 + 0] * aux_g[1 * objects_count * 3 + j * 3 + 2]
                + aux_c[2 * dim_nodes_minus_1 + 0] * aux_g[2 * objects_count * 3 + j * 3 + 2]
                + aux_c[3 * dim_nodes_minus_1 + 0] * aux_g[3 * objects_count * 3 + j * 3 + 2]
                + aux_c[4 * dim_nodes_minus_1 + 0] * aux_g[4 * objects_count * 3 + j * 3 + 2]
                + aux_c[5 * dim_nodes_minus_1 + 0] * aux_g[5 * objects_count * 3 + j * 3 + 2]
                + aux_c[6 * dim_nodes_minus_1 + 0] * aux_g[6 * objects_count * 3 + j * 3 + 2]
            );
        }
        else
        {
            continue;
        }

        if (i >= 2) 
        {
            aux_b[1 * objects_count * 3 + j * 3 + 0] = (
                aux_c[1 * dim_nodes_minus_1 + 1] * aux_g[1 * objects_count * 3 + j * 3 + 0]
                + aux_c[2 * dim_nodes_minus_1 + 1] * aux_g[2 * objects_count * 3 + j * 3 + 0]
                + aux_c[3 * dim_nodes_minus_1 + 1] * aux_g[3 * objects_count * 3 + j * 3 + 0]
                + aux_c[4 * dim_nodes_minus_1 + 1] * aux_g[4 * objects_count * 3 + j * 3 + 0]
                + aux_c[5 * dim_nodes_minus_1 + 1] * aux_g[5 * objects_count * 3 + j * 3 + 0]
                + aux_c[6 * dim_nodes_minus_1 + 1] * aux_g[6 * objects_count * 3 + j * 3 + 0]
            );
            aux_b[1 * objects_count * 3 + j * 3 + 1] = (
                aux_c[1 * dim_nodes_minus_1 + 1] * aux_g[1 * objects_count * 3 + j * 3 + 1]
                + aux_c[2 * dim_nodes_minus_1 + 1] * aux_g[2 * objects_count * 3 + j * 3 + 1]
                + aux_c[3 * dim_nodes_minus_1 + 1] * aux_g[3 * objects_count * 3 + j * 3 + 1]
                + aux_c[4 * dim_nodes_minus_1 + 1] * aux_g[4 * objects_count * 3 + j * 3 + 1]
                + aux_c[5 * dim_nodes_minus_1 + 1] * aux_g[5 * objects_count * 3 + j * 3 + 1]
                + aux_c[6 * dim_nodes_minus_1 + 1] * aux_g[6 * objects_count * 3 + j * 3 + 1]
            );
            aux_b[1 * objects_count * 3 + j * 3 + 2] = (
                aux_c[1 * dim_nodes_minus_1 + 1] * aux_g[1 * objects_count * 3 + j * 3 + 2]
                + aux_c[2 * dim_nodes_minus_1 + 1] * aux_g[2 * objects_count * 3 + j * 3 + 2]
                + aux_c[3 * dim_nodes_minus_1 + 1] * aux_g[3 * objects_count * 3 + j * 3 + 2]
                + aux_c[4 * dim_nodes_minus_1 + 1] * aux_g[4 * objects_count * 3 + j * 3 + 2]
                + aux_c[5 * dim_nodes_minus_1 + 1] * aux_g[5 * objects_count * 3 + j * 3 + 2]
                + aux_c[6 * dim_nodes_minus_1 + 1] * aux_g[6 * objects_count * 3 + j * 3 + 2]
            );
        }
        else
        {
            continue;
        }

        if (i >= 3) 
        {
            aux_b[2 * objects_count * 3 + j * 3 + 0] = (
                aux_c[2 * dim_nodes_minus_1 + 2] * aux_g[2 * objects_count * 3 + j * 3 + 0]
                + aux_c[3 * dim_nodes_minus_1 + 2] * aux_g[3 * objects_count * 3 + j * 3 + 0]
                + aux_c[4 * dim_nodes_minus_1 + 2] * aux_g[4 * objects_count * 3 + j * 3 + 0]
                + aux_c[5 * dim_nodes_minus_1 + 2] * aux_g[5 * objects_count * 3 + j * 3 + 0]
                + aux_c[6 * dim_nodes_minus_1 + 2] * aux_g[6 * objects_count * 3 + j * 3 + 0]
            );
            aux_b[2 * objects_count * 3 + j * 3 + 1] = (
                aux_c[2 * dim_nodes_minus_1 + 2] * aux_g[2 * objects_count * 3 + j * 3 + 1]
                + aux_c[3 * dim_nodes_minus_1 + 2] * aux_g[3 * objects_count * 3 + j * 3 + 1]
                + aux_c[4 * dim_nodes_minus_1 + 2] * aux_g[4 * objects_count * 3 + j * 3 + 1]
                + aux_c[5 * dim_nodes_minus_1 + 2] * aux_g[5 * objects_count * 3 + j * 3 + 1]
                + aux_c[6 * dim_nodes_minus_1 + 2] * aux_g[6 * objects_count * 3 + j * 3 + 1]
            );
            aux_b[2 * objects_count * 3 + j * 3 + 2] = (
                aux_c[2 * dim_nodes_minus_1 + 2] * aux_g[2 * objects_count * 3 + j * 3 + 2]
                + aux_c[3 * dim_nodes_minus_1 + 2] * aux_g[3 * objects_count * 3 + j * 3 + 2]
                + aux_c[4 * dim_nodes_minus_1 + 2] * aux_g[4 * objects_count * 3 + j * 3 + 2]
                + aux_c[5 * dim_nodes_minus_1 + 2] * aux_g[5 * objects_count * 3 + j * 3 + 2]
                + aux_c[6 * dim_nodes_minus_1 + 2] * aux_g[6 * objects_count * 3 + j * 3 + 2]
            );
        }
        else
        {
            continue;
        }

        if (i >= 4) 
        {
            aux_b[3 * objects_count * 3 + j * 3 + 0] = (
                aux_c[3 * dim_nodes_minus_1 + 3] * aux_g[3 * objects_count * 3 + j * 3 + 0]
                + aux_c[4 * dim_nodes_minus_1 + 3] * aux_g[4 * objects_count * 3 + j * 3 + 0]
                + aux_c[5 * dim_nodes_minus_1 + 3] * aux_g[5 * objects_count * 3 + j * 3 + 0]
                + aux_c[6 * dim_nodes_minus_1 + 3] * aux_g[6 * objects_count * 3 + j * 3 + 0]
            );
            aux_b[3 * objects_count * 3 + j * 3 + 1] = (
                aux_c[3 * dim_nodes_minus_1 + 3] * aux_g[3 * objects_count * 3 + j * 3 + 1]
                + aux_c[4 * dim_nodes_minus_1 + 3] * aux_g[4 * objects_count * 3 + j * 3 + 1]
                + aux_c[5 * dim_nodes_minus_1 + 3] * aux_g[5 * objects_count * 3 + j * 3 + 1]
                + aux_c[6 * dim_nodes_minus_1 + 3] * aux_g[6 * objects_count * 3 + j * 3 + 1]
            );
            aux_b[3 * objects_count * 3 + j * 3 + 2] = (
                aux_c[3 * dim_nodes_minus_1 + 3] * aux_g[3 * objects_count * 3 + j * 3 + 2]
                + aux_c[4 * dim_nodes_minus_1 + 3] * aux_g[4 * objects_count * 3 + j * 3 + 2]
                + aux_c[5 * dim_nodes_minus_1 + 3] * aux_g[5 * objects_count * 3 + j * 3 + 2]
                + aux_c[6 * dim_nodes_minus_1 + 3] * aux_g[6 * objects_count * 3 + j * 3 + 2]
            );
        }
        else
        {
            continue;
        }

        if (i >= 5)
        {
            aux_b[4 * objects_count * 3 + j * 3 + 0] = (
                aux_c[4 * dim_nodes_minus_1 + 4] * aux_g[4 * objects_count * 3 + j * 3 + 0]
                + aux_c[5 * dim_nodes_minus_1 + 4] * aux_g[5 * objects_count * 3 + j * 3 + 0]
                + aux_c[6 * dim_nodes_minus_1 + 4] * aux_g[6 * objects_count * 3 + j * 3 + 0]
            );
            aux_b[4 * objects_count * 3 + j * 3 + 1] = (
                aux_c[4 * dim_nodes_minus_1 + 4] * aux_g[4 * objects_count * 3 + j * 3 + 1]
                + aux_c[5 * dim_nodes_minus_1 + 4] * aux_g[5 * objects_count * 3 + j * 3 + 1]
                + aux_c[6 * dim_nodes_minus_1 + 4] * aux_g[6 * objects_count * 3 + j * 3 + 1]
            );
            aux_b[4 * objects_count * 3 + j * 3 + 2] = (
                aux_c[4 * dim_nodes_minus_1 + 4] * aux_g[4 * objects_count * 3 + j * 3 + 2]
                + aux_c[5 * dim_nodes_minus_1 + 4] * aux_g[5 * objects_count * 3 + j * 3 + 2]
                + aux_c[6 * dim_nodes_minus_1 + 4] * aux_g[6 * objects_count * 3 + j * 3 + 2]
            );
        }
        else
        {
            continue;
        }

        if (i >= 6)
        {
            aux_b[5 * objects_count * 3 + j * 3 + 0] = (
                aux_c[5 * dim_nodes_minus_1 + 5] * aux_g[5 * objects_count * 3 + j * 3 + 0]
                + aux_c[6 * dim_nodes_minus_1 + 5] * aux_g[6 * objects_count * 3 + j * 3 + 0]
            );
            aux_b[5 * objects_count * 3 + j * 3 + 1] = (
                aux_c[5 * dim_nodes_minus_1 + 5] * aux_g[5 * objects_count * 3 + j * 3 + 1]
                + aux_c[6 * dim_nodes_minus_1 + 5] * aux_g[6 * objects_count * 3 + j * 3 + 1]
            );
            aux_b[5 * objects_count * 3 + j * 3 + 2] = (
                aux_c[5 * dim_nodes_minus_1 + 5] * aux_g[5 * objects_count * 3 + j * 3 + 2]
                + aux_c[6 * dim_nodes_minus_1 + 5] * aux_g[6 * objects_count * 3 + j * 3 + 2]
            );
        }
        else
        {
            continue;
        }

        if (i >= 7)
        {
            aux_b[6 * objects_count * 3 + j * 3 + 0] = (
                aux_c[6 * dim_nodes_minus_1 + 6] * aux_g[6 * objects_count * 3 + j * 3 + 0]
            );
            aux_b[6 * objects_count * 3 + j * 3 + 1] = (
                aux_c[6 * dim_nodes_minus_1 + 6] * aux_g[6 * objects_count * 3 + j * 3 + 1]
            );
            aux_b[6 * objects_count * 3 + j * 3 + 2] = (
                aux_c[6 * dim_nodes_minus_1 + 6] * aux_g[6 * objects_count * 3 + j * 3 + 2]
            );
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
    memcpy(F, aux_a, (i + 1) * objects_count * 3 * sizeof(real));
    
    // Update aux_g
    for (int j = 0; j < objects_count; j++)
    {
        if (i >= 1)
        {
            aux_g[0 * objects_count * 3 + j * 3 + 0] = (
                (F[1 * objects_count * 3 + j * 3 + 0] - F[0 * objects_count * 3 + j * 3 + 0]) * aux_r[1 * dim_nodes + 0]
            );

            aux_g[0 * objects_count * 3 + j * 3 + 1] = (
                (F[1 * objects_count * 3 + j * 3 + 1] - F[0 * objects_count * 3 + j * 3 + 1]) * aux_r[1 * dim_nodes + 0]
            );
            
            aux_g[0 * objects_count * 3 + j * 3 + 2] = (
                (F[1 * objects_count * 3 + j * 3 + 2] - F[0 * objects_count * 3 + j * 3 + 2]) * aux_r[1 * dim_nodes + 0]
            );
        }
        else
        {
            continue;
        }
            
        if (i >= 2)
        {
            aux_g[1 * objects_count * 3 + j * 3 + 0] = (
                ((F[2 * objects_count * 3 + j * 3 + 0] - F[0 * objects_count * 3 + j * 3 + 0]) * aux_r[2 * dim_nodes + 0] 
                - aux_g[0 * objects_count * 3 + j * 3 + 0]) * aux_r[2 * dim_nodes + 1]
            );

            aux_g[1 * objects_count * 3 + j * 3 + 1] = (
                ((F[2 * objects_count * 3 + j * 3 + 1] - F[0 * objects_count * 3 + j * 3 + 1]) * aux_r[2 * dim_nodes + 0] 
                - aux_g[0 * objects_count * 3 + j * 3 + 1]) * aux_r[2 * dim_nodes + 1]
            );

            aux_g[1 * objects_count * 3 + j * 3 + 2] = (
                ((F[2 * objects_count * 3 + j * 3 + 2] - F[0 * objects_count * 3 + j * 3 + 2]) * aux_r[2 * dim_nodes + 0] 
                - aux_g[0 * objects_count * 3 + j * 3 + 2]) * aux_r[2 * dim_nodes + 1]
            );
        }
        else 
        {
            continue;
        }

        if (i >= 3)
        {
            aux_g[2 * objects_count * 3 + j * 3 + 0] = (
                ((F[3 * objects_count * 3 + j * 3 + 0] - F[0 * objects_count * 3 + j * 3 + 0]) * aux_r[3 * dim_nodes + 0] 
                - aux_g[0 * objects_count * 3 + j * 3 + 0]) * aux_r[3 * dim_nodes + 1] 
                - aux_g[1 * objects_count * 3 + j * 3 + 0]
            ) * aux_r[3 * dim_nodes + 2];

            aux_g[2 * objects_count * 3 + j * 3 + 1] = (
                ((F[3 * objects_count * 3 + j * 3 + 1] - F[0 * objects_count * 3 + j * 3 + 1]) * aux_r[3 * dim_nodes + 0] 
                - aux_g[0 * objects_count * 3 + j * 3 + 1]) * aux_r[3 * dim_nodes + 1] 
                - aux_g[1 * objects_count * 3 + j * 3 + 1]
            ) * aux_r[3 * dim_nodes + 2];

            aux_g[2 * objects_count * 3 + j * 3 + 2] = (
                ((F[3 * objects_count * 3 + j * 3 + 2] - F[0 * objects_count * 3 + j * 3 + 2]) * aux_r[3 * dim_nodes + 0] 
                - aux_g[0 * objects_count * 3 + j * 3 + 2]) * aux_r[3 * dim_nodes + 1] 
                - aux_g[1 * objects_count * 3 + j * 3 + 2]
            ) * aux_r[3 * dim_nodes + 2];
        }
        else
        {
            continue;
        }
            
        if (i >= 4)
        {
            aux_g[3 * objects_count * 3 + j * 3 + 0] = (
                (((F[4 * objects_count * 3 + j * 3 + 0] - F[0 * objects_count * 3 + j * 3 + 0]) * aux_r[4 * dim_nodes + 0] 
                - aux_g[0 * objects_count * 3 + j * 3 + 0]) * aux_r[4 * dim_nodes + 1] - aux_g[1 * objects_count * 3 + j * 3 + 0])
                * aux_r[4 * dim_nodes + 2]
                - aux_g[2 * objects_count * 3 + j * 3 + 0]
            ) * aux_r[4 * dim_nodes + 3];

            aux_g[3 * objects_count * 3 + j * 3 + 1] = (
                (((F[4 * objects_count * 3 + j * 3 + 1] - F[0 * objects_count * 3 + j * 3 + 1]) * aux_r[4 * dim_nodes + 0] 
                - aux_g[0 * objects_count * 3 + j * 3 + 1]) * aux_r[4 * dim_nodes + 1] - aux_g[1 * objects_count * 3 + j * 3 + 1])
                * aux_r[4 * dim_nodes + 2]
                - aux_g[2 * objects_count * 3 + j * 3 + 1]
            ) * aux_r[4 * dim_nodes + 3];

            aux_g[3 * objects_count * 3 + j * 3 + 2] = (
                (((F[4 * objects_count * 3 + j * 3 + 2] - F[0 * objects_count * 3 + j * 3 + 2]) * aux_r[4 * dim_nodes + 0] 
                - aux_g[0 * objects_count * 3 + j * 3 + 2]) * aux_r[4 * dim_nodes + 1] - aux_g[1 * objects_count * 3 + j * 3 + 2])
                * aux_r[4 * dim_nodes + 2]
                - aux_g[2 * objects_count * 3 + j * 3 + 2]
            ) * aux_r[4 * dim_nodes + 3];
        }
        else
        {
            continue;
        }

        if (i >= 5)
        {
            aux_g[4 * objects_count * 3 + j * 3 + 0] = (
                (
                    (((F[5 * objects_count * 3 + j * 3 + 0] - F[0 * objects_count * 3 + j * 3 + 0]) * aux_r[5 * dim_nodes + 0] 
                    - aux_g[0 * objects_count * 3 + j * 3 + 0]) * aux_r[5 * dim_nodes + 1] 
                    - aux_g[1 * objects_count * 3 + j * 3 + 0])
                    * aux_r[5 * dim_nodes + 2]
                    - aux_g[2 * objects_count * 3 + j * 3 + 0]
                )
                * aux_r[5 * dim_nodes + 3]
                - aux_g[3 * objects_count * 3 + j * 3 + 0]
            ) * aux_r[5 * dim_nodes + 4];

            aux_g[4 * objects_count * 3 + j * 3 + 1] = (
                (
                    (((F[5 * objects_count * 3 + j * 3 + 1] - F[0 * objects_count * 3 + j * 3 + 1]) * aux_r[5 * dim_nodes + 0] 
                    - aux_g[0 * objects_count * 3 + j * 3 + 1]) * aux_r[5 * dim_nodes + 1] 
                    - aux_g[1 * objects_count * 3 + j * 3 + 1])
                    * aux_r[5 * dim_nodes + 2]
                    - aux_g[2 * objects_count * 3 + j * 3 + 1]
                )
                * aux_r[5 * dim_nodes + 3]
                - aux_g[3 * objects_count * 3 + j * 3 + 1]
            ) * aux_r[5 * dim_nodes + 4];

            aux_g[4 * objects_count * 3 + j * 3 + 2] = (
                (
                    (((F[5 * objects_count * 3 + j * 3 + 2] - F[0 * objects_count * 3 + j * 3 + 2]) * aux_r[5 * dim_nodes + 0] 
                    - aux_g[0 * objects_count * 3 + j * 3 + 2]) * aux_r[5 * dim_nodes + 1] 
                    - aux_g[1 * objects_count * 3 + j * 3 + 2])
                    * aux_r[5 * dim_nodes + 2]
                    - aux_g[2 * objects_count * 3 + j * 3 + 2]
                )
                * aux_r[5 * dim_nodes + 3]
                - aux_g[3 * objects_count * 3 + j * 3 + 2]
            ) * aux_r[5 * dim_nodes + 4];
        }
        else
        {
            continue;
        }

        if (i >= 6)
        {
            aux_g[5 * objects_count * 3 + j * 3 + 0] = (
                (
                    (
                        (((F[6 * objects_count * 3 + j * 3 + 0] - F[0 * objects_count * 3 + j * 3 + 0]) * aux_r[6 * dim_nodes + 0] 
                        - aux_g[0 * objects_count * 3 + j * 3 + 0]) * aux_r[6 * dim_nodes + 1] 
                        - aux_g[1 * objects_count * 3 + j * 3 + 0])
                        * aux_r[6 * dim_nodes + 2]
                        - aux_g[2 * objects_count * 3 + j * 3 + 0]
                    )
                    * aux_r[6 * dim_nodes + 3]
                    - aux_g[3 * objects_count * 3 + j * 3 + 0]
                )
                * aux_r[6 * dim_nodes + 4]
                - aux_g[4 * objects_count * 3 + j * 3 + 0]
            ) * aux_r[6 * dim_nodes + 5];

            aux_g[5 * objects_count * 3 + j * 3 + 1] = (
                (
                    (
                        (((F[6 * objects_count * 3 + j * 3 + 1] - F[0 * objects_count * 3 + j * 3 + 1]) * aux_r[6 * dim_nodes + 0] 
                        - aux_g[0 * objects_count * 3 + j * 3 + 1]) * aux_r[6 * dim_nodes + 1] 
                        - aux_g[1 * objects_count * 3 + j * 3 + 1])
                        * aux_r[6 * dim_nodes + 2]
                        - aux_g[2 * objects_count * 3 + j * 3 + 1]
                    )
                    * aux_r[6 * dim_nodes + 3]
                    - aux_g[3 * objects_count * 3 + j * 3 + 1]
                )
                * aux_r[6 * dim_nodes + 4]
                - aux_g[4 * objects_count * 3 + j * 3 + 1]
            ) * aux_r[6 * dim_nodes + 5];

            aux_g[5 * objects_count * 3 + j * 3 + 2] = (
                (
                    (
                        (((F[6 * objects_count * 3 + j * 3 + 2] - F[0 * objects_count * 3 + j * 3 + 2]) * aux_r[6 * dim_nodes + 0] 
                        - aux_g[0 * objects_count * 3 + j * 3 + 2]) * aux_r[6 * dim_nodes + 1] 
                        - aux_g[1 * objects_count * 3 + j * 3 + 2])
                        * aux_r[6 * dim_nodes + 2]
                        - aux_g[2 * objects_count * 3 + j * 3 + 2]
                    )
                    * aux_r[6 * dim_nodes + 3]
                    - aux_g[3 * objects_count * 3 + j * 3 + 2]
                )
                * aux_r[6 * dim_nodes + 4]
                - aux_g[4 * objects_count * 3 + j * 3 + 2]
            ) * aux_r[6 * dim_nodes + 5];
        }
        else
        {
            continue;
        }

        if (i >= 7)
        {
            aux_g[6 * objects_count * 3 + j * 3 + 0] = (
                (
                    (
                        (
                            (((F[7 * objects_count * 3 + j * 3 + 0] - F[0 * objects_count * 3 + j * 3 + 0]) * aux_r[7 * dim_nodes + 0] - aux_g[0 * objects_count * 3 + j * 3 + 0]) 
                            * aux_r[7 * dim_nodes + 1] 
                            - aux_g[1 * objects_count * 3 + j * 3 + 0])
                            * aux_r[7 * dim_nodes + 2]
                            - aux_g[2 * objects_count * 3 + j * 3 + 0]
                        )
                        * aux_r[7 * dim_nodes + 3]
                        - aux_g[3 * objects_count * 3 + j * 3 + 0]
                    )
                    * aux_r[7 * dim_nodes + 4]
                    - aux_g[4 * objects_count * 3 + j * 3 + 0]
                )
                * aux_r[7 * dim_nodes + 5]
                - aux_g[5 * objects_count * 3 + j * 3 + 0]
            ) * aux_r[7 * dim_nodes + 6];          

            aux_g[6 * objects_count * 3 + j * 3 + 1] = (
                (
                    (
                        (
                            (((F[7 * objects_count * 3 + j * 3 + 1] - F[0 * objects_count * 3 + j * 3 + 1]) * aux_r[7 * dim_nodes + 0] - aux_g[0 * objects_count * 3 + j * 3 + 1]) 
                            * aux_r[7 * dim_nodes + 1] 
                            - aux_g[1 * objects_count * 3 + j * 3 + 1])
                            * aux_r[7 * dim_nodes + 2]
                            - aux_g[2 * objects_count * 3 + j * 3 + 1]
                        )
                        * aux_r[7 * dim_nodes + 3]
                        - aux_g[3 * objects_count * 3 + j * 3 + 1]
                    )
                    * aux_r[7 * dim_nodes + 4]
                    - aux_g[4 * objects_count * 3 + j * 3 + 1]
                )
                * aux_r[7 * dim_nodes + 5]
                - aux_g[5 * objects_count * 3 + j * 3 + 1]
            ) * aux_r[7 * dim_nodes + 6];       

            aux_g[6 * objects_count * 3 + j * 3 + 2] = (
                (
                    (
                        (
                            (((F[7 * objects_count * 3 + j * 3 + 2] - F[0 * objects_count * 3 + j * 3 + 2]) * aux_r[7 * dim_nodes + 0] - aux_g[0 * objects_count * 3 + j * 3 + 2]) 
                            * aux_r[7 * dim_nodes + 1] 
                            - aux_g[1 * objects_count * 3 + j * 3 + 2])
                            * aux_r[7 * dim_nodes + 2]
                            - aux_g[2 * objects_count * 3 + j * 3 + 2]
                        )
                        * aux_r[7 * dim_nodes + 3]
                        - aux_g[3 * objects_count * 3 + j * 3 + 2]
                    )
                    * aux_r[7 * dim_nodes + 4]
                    - aux_g[4 * objects_count * 3 + j * 3 + 2]
                )
                * aux_r[7 * dim_nodes + 5]
                - aux_g[5 * objects_count * 3 + j * 3 + 2]
            ) * aux_r[7 * dim_nodes + 6];             
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
                for (int k = 0; k < 3; k++)
                {
                    delta_aux_b[i * objects_count * 3 + j * 3 + k] = (
                        aux_b[i * objects_count * 3 + j * 3 + k] - aux_e[i * objects_count * 3 + j * 3 + k]
                    );
                }
            }
        }
    }
    else
    {
        // Empty delta_aux_b
        for (int i = 0; i < dim_nodes_minus_1 * objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                delta_aux_b[i * 3 + j] = 0.0;
            }
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
            for (int k = 0; k < 3; k++)
            {
                aux_b[i * objects_count * 3 + j * 3 + k] = (
                    aux_e[i * objects_count * 3 + j * 3 + k] + delta_aux_b[i * objects_count * 3 + j * 3 + k]
                );
            }
        }
    }
}
