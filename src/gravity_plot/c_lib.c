#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#ifdef WIN32DLL_EXPORTS
    #define WIN32DLL_API __declspec(dllexport)
#else
    #define WIN32DLL_API 
#endif

#define NPTS 50000

typedef int64_t int64;
typedef double real;

typedef struct Solutions 
{
    double* sol_state;
    double* sol_time;
    double* sol_dt;
    double* m;
    double G;
    int objects_count;
} Solutions;

void free_memory_real(real *ptr);
void initialize_system(
    const char *restrict system,
    real **x,
    real **v,
    real **m,
    int *restrict objects_count,
    real *restrict G,
    int *restrict custom_sys_flag
);
real abs_max_vec(const real *restrict vec, int vec_length);
real vec_norm(const real *restrict vec, int vec_length);
void compute_energy(
    int objects_count, 
    int npts,
    int *restrict count, 
    double *restrict energy, 
    const double (*restrict sol_state)[objects_count * 6], 
    const double *restrict m, 
    real G
);
void acceleration(int objects_count, real *restrict x, real *restrict a, const real *restrict m, real G);
Solutions euler(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
);
Solutions euler_cromer(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
);
Solutions rk4(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
);
Solutions leapfrog(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
);
Solutions rk_embedded(
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
    int custom_sys_objects_count
);
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
void rk_embedded_butcher_tableaus(
    int order,
    int *restrict power,
    int *restrict power_test,
    real **coeff,
    int *restrict len_weights,
    real **weights,
    real **weights_test
);
Solutions ias15(
    const char *restrict system,
    double *restrict t,
    double tf, 
    double input_tolerance,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
);
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
    real *restrict temp_v_err_comp_sum
);
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
void ias15_approx_vel_step(
    int objects_count,
    real *restrict v,
    real *restrict v0,
    real *restrict a0,
    real *restrict aux_b,
    real dt,
    real *restrict temp_v_err_comp_sum
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
void ias15_radau_spacing(real *restrict nodes);
void ias15_aux_c(real *restrict aux_c);
void ias15_aux_r(real *aux_r);
real ias15_initial_dt(
    int power,
    real *restrict x,
    real *restrict v,
    real *restrict a,
    const real *m,
    int objects_count,
    real G
);

// Free memory in type real
WIN32DLL_API void free_memory_real(real *ptr)
{
    free(ptr);
}

WIN32DLL_API void initialize_system(
    const char *restrict system,
    real **x,
    real **v,
    real **m,
    int *restrict objects_count,
    real *restrict G,
    int *restrict custom_sys_flag
)
{
    // Conversion factor from km^3 s^-2 to AU^3 d^-2
    real CONVERSION_FACTOR = ((real) 86400.0L * 86400.0L) / (149597870.7L * 149597870.7L * 149597870.7L);
    // GM values (km^3 s^-2)
    // ref: https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf
    const real GM_SI_SUN = 132712440041.279419L;
    const real GM_SI_MERCURY = 22031.868551L;
    const real GM_SI_VENUS = 324858.592000L;
    const real GM_SI_EARTH = 398600.435507L;
    const real GM_SI_MARS = 42828.375816L;
    const real GM_SI_JUPITER = 126712764.100000L;
    const real GM_SI_SATURN = 37940584.841800L;
    const real GM_SI_URANUS = 5794556.400000L;
    const real GM_SI_NEPTUNE = 6836527.100580L;
    const real GM_SI_MOON = 4902.800118L;
    const real GM_SI_PLUTO = 975.500000L;
    const real GM_SI_CERES = 62.62890L;
    const real GM_SI_VESTA = 17.288245L;

    // GM values (AU^3 d^-2)
    const real GM_SUN = 132712440041.279419L * CONVERSION_FACTOR;
    // const real GM_MERCURY = 22031.868551L * CONVERSION_FACTOR;
    // const real GM_VENUS = 324858.592000L * CONVERSION_FACTOR;
    // const real GM_EARTH = 398600.435507L * CONVERSION_FACTOR;
    // const real GM_MARS = 42828.375816L * CONVERSION_FACTOR;
    // const real GM_JUPITER = 126712764.100000L * CONVERSION_FACTOR;
    // const real GM_SATURN = 37940584.841800L * CONVERSION_FACTOR;
    // const real GM_URANUS = 5794556.400000L * CONVERSION_FACTOR;
    // const real GM_NEPTUNE = 6836527.100580L * CONVERSION_FACTOR;
    // const real GM_MOON = 4902.800118L * CONVERSION_FACTOR;
    // const real GM_PLUTO = 975.500000L * CONVERSION_FACTOR;
    // const real GM_CERES = 62.62890L * CONVERSION_FACTOR;
    // const real GM_VESTA = 17.288245L * CONVERSION_FACTOR;

    // Solar system masses (M_sun^-1)
    const real MASS_SUN = 1.0;
    const real MASS_MERCURY = GM_SI_MERCURY / GM_SI_SUN;
    const real MASS_VENUS = GM_SI_VENUS / GM_SI_SUN;
    const real MASS_EARTH = GM_SI_EARTH / GM_SI_SUN;
    const real MASS_MARS = GM_SI_MARS / GM_SI_SUN;
    const real MASS_JUPITER = GM_SI_JUPITER / GM_SI_SUN;
    const real MASS_SATURN = GM_SI_SATURN / GM_SI_SUN;
    const real MASS_URANUS = GM_SI_URANUS / GM_SI_SUN;
    const real MASS_NEPTUNE = GM_SI_NEPTUNE / GM_SI_SUN;
    const real MASS_MOON = GM_SI_MOON / GM_SI_SUN;
    const real MASS_PLUTO = GM_SI_PLUTO / GM_SI_SUN;
    const real MASS_CERES = GM_SI_CERES / GM_SI_SUN;
    const real MASS_VESTA = GM_SI_VESTA / GM_SI_SUN;

    // Gravitational constant (kg^-1 m^3 s^-2):
    // const real G_SI = 6.67430e-11;
    // Gravitational constant (M_sun^-1 AU^3 d^-2):
    *G = GM_SUN;

    /*
    * Solar system position and velocities data
    * Units: AU-D
    * Coordinate center: Solar System Barycenter
    * Data dated on A.D. 2024-Jan-01 00:00:00.0000 TDB
    * Computational data generated by NASA JPL Horizons System https://ssd.jpl.nasa.gov/horizons/
    */
    const real POS_SUN[3] = {-7.967955691533730e-03L, -2.906227441573178e-03L, 2.103054301547123e-04L};
    const real POS_MERCURY[3] = {-2.825983269538632e-01L, 1.974559795958082e-01L, 4.177433558063677e-02L};
    const real POS_VENUS[3] = {-7.232103701666379e-01L, -7.948302026312400e-02L, 4.042871428174315e-02L};
    const real POS_EARTH[3] = {-1.738192017257054e-01L, 9.663245550235138e-01L, 1.553901854897183e-04L};
    const real POS_MARS[3] = {-3.013262392582653e-01L, -1.454029331393295e00L, -2.300531433991428e-02L};
    const real POS_JUPITER[3] = {3.485202469657674e00L, 3.552136904413157e00L, -9.271035442798399e-02L};
    const real POS_SATURN[3] = {8.988104223143450e00L, -3.719064854634689e00L, -2.931937777323593e-01L};
    const real POS_URANUS[3] = {1.226302417897505e01L, 1.529738792480545e01L, -1.020549026883563e-01L};
    const real POS_NEPTUNE[3] = {2.983501460984741e01L, -1.793812957956852e00L, -6.506401132254588e-01L};
    const real POS_MOON[3] = {-1.762788124769829e-01L, 9.674377513177153e-01L, 3.236901585768862e-04L};
    const real POS_PLUTO[3] = {1.720200478843485e01L, -3.034155683573043e01L, -1.729127607100611e00L};
    const real POS_CERES[3] = {-1.103880510367569e00L, -2.533340440444230e00L, 1.220283937721780e-01L};
    const real POS_VESTA[3] = {-8.092549658731499e-02L, 2.558381434460076e00L, -6.695836142398572e-02L};

    const real VEL_SUN[3] = {4.875094764261564e-06L, -7.057133213976680e-06L, -4.573453713094512e-08L};
    const real VEL_MERCURY[3] = {-2.232165900189702e-02L, -2.157207103176252e-02L, 2.855193410495743e-04L};
    const real VEL_VENUS[3] = {2.034068201002341e-03L, -2.020828626592994e-02L, -3.945639843855159e-04L};
    const real VEL_EARTH[3] = {-1.723001232538228e-02L, -2.967721342618870e-03L, 6.382125383116755e-07L};
    const real VEL_MARS[3] = {1.424832259345280e-02L, -1.579236181580905e-03L, -3.823722796161561e-04L};
    const real VEL_JUPITER[3] = {-5.470970658852281e-03L, 5.642487338479145e-03L, 9.896190602066252e-05L};
    const real VEL_SATURN[3] = {1.822013845554067e-03L, 5.143470425888054e-03L, -1.617235904887937e-04L};
    const real VEL_URANUS[3] = {-3.097615358317413e-03L, 2.276781932345769e-03L, 4.860433222241686e-05L};
    const real VEL_NEPTUNE[3] = {1.676536611817232e-04L, 3.152098732861913e-03L, -6.877501095688201e-05L};
    const real VEL_MOON[3] = {-1.746667306153906e-02L, -3.473438277358121e-03L, -3.359028758606074e-05L};
    const real VEL_PLUTO[3] = {2.802810313667557e-03L, 8.492056438614633e-04L, -9.060790113327894e-04L};
    const real VEL_CERES[3] = {8.978653480111301e-03L, -4.873256528198994e-03L, -1.807162046049230e-03L};
    const real VEL_VESTA[3] = {-1.017876585480054e-02L, -5.452367109338154e-04L, 1.255870551153315e-03L};
    

    // Pre-defined systems
    if (strcmp(system, "circular_binary_orbit") == 0) 
    {
        *objects_count = 2;
        *x = malloc(*objects_count * 3 * sizeof(real));
        *v = malloc(*objects_count * 3 * sizeof(real));
        *m = malloc(*objects_count * sizeof(real));

        (*x)[0] = 1.0;
        (*x)[1] = 0.0;
        (*x)[2] = 0.0;

        (*x)[3] = -1.0;
        (*x)[4] = 0.0;
        (*x)[5] = 0.0; 

        (*v)[0] = 0.0;
        (*v)[1] = 0.5;
        (*v)[2] = 0.0;

        (*v)[3] = 0.0;
        (*v)[4] = -0.5;
        (*v)[5] = 0.0;
        
        (*m)[0] = 1.0 / *G;
        (*m)[1] = 1.0 / *G;
    }
    else if (strcmp(system, "eccentric_binary_orbit") == 0) 
    {
        *objects_count = 2;
        *x = malloc(*objects_count * 3 * sizeof(real));
        *v = malloc(*objects_count * 3 * sizeof(real));
        *m = malloc(*objects_count * sizeof(real));

        (*x)[0] = 1.0;
        (*x)[1] = 0.0;
        (*x)[2] = 0.0;

        (*x)[3] = -1.25;
        (*x)[4] = 0.0;
        (*x)[5] = 0.0; 

        (*v)[0] = 0.0;
        (*v)[1] = 0.5;
        (*v)[2] = 0.0;

        (*v)[3] = 0.0;
        (*v)[4] = -0.625;
        (*v)[5] = 0.0;
        
        (*m)[0] = 1.0 / *G;
        (*m)[1] = 0.8 / *G;
    }
    else if (strcmp(system, "3d_helix") == 0) 
    {
        *objects_count = 3;
        *x = malloc(*objects_count * 3 * sizeof(real));
        *v = malloc(*objects_count * 3 * sizeof(real));
        *m = malloc(*objects_count * sizeof(real));
        
        (*x)[0] = 0.0;
        (*x)[1] = 0.0;
        (*x)[2] = -1.0;

        (*x)[3] = -sqrt(3.0) / 2.0;
        (*x)[4] = 0.0;
        (*x)[5] = 0.5;

        (*x)[6] = sqrt(3.0) / 2.0;
        (*x)[7] = 0.0;
        (*x)[8] = 0.5;

        real v0 = sqrt(1.0 / sqrt(3));

        (*v)[0] = -v0;
        (*v)[1] = 0.5;
        (*v)[2] = 0.0;

        (*v)[3] = 0.5 * v0;
        (*v)[4] = 0.5;
        (*v)[5] = (sqrt(3.0) / 2.0) * v0;

        (*v)[6] = 0.5 * v0;
        (*v)[7] = 0.5;
        (*v)[8] = -(sqrt(3.0) / 2.0) * v0;
        
        (*m)[0] = 1.0 / *G;
        (*m)[1] = 1.0 / *G;
        (*m)[2] = 1.0 / *G;
    }
    else if (strcmp(system, "sun_earth_moon") == 0) 
    {
        *objects_count = 3;
        *x = malloc(*objects_count * 3 * sizeof(real));
        *v = malloc(*objects_count * 3 * sizeof(real));
        *m = malloc(*objects_count * sizeof(real));

        (*m)[0] = MASS_SUN;
        (*m)[1] = MASS_EARTH;
        (*m)[2] = MASS_MOON;
        
        real R_CM[3];
        real V_CM[3];
        const real M = (*m)[0] + (*m)[1] + (*m)[2];
        for (int i = 0; i < 3; i++)
        {
            R_CM[i] = 1 / M * ((*m)[0] * POS_SUN[i] + (*m)[1] * POS_EARTH[i] + (*m)[2] * POS_MOON[i]);
            V_CM[i] = 1 / M * ((*m)[0] * VEL_SUN[i] + (*m)[1] * VEL_EARTH[i] + (*m)[2] * VEL_MOON[i]);
        }

        for (int i = 0; i < 3; i++)
        {
            (*x)[0 * 3 + i] = POS_SUN[i] - R_CM[i];
            (*x)[1 * 3 + i] = POS_EARTH[i] - R_CM[i];
            (*x)[2 * 3 + i] = POS_MOON[i] - R_CM[i];

            (*v)[0 * 3 + i] = VEL_SUN[i] - V_CM[i];
            (*v)[1 * 3 + i] = VEL_EARTH[i] - V_CM[i];
            (*v)[2 * 3 + i] = VEL_MOON[i] - V_CM[i];
        }
    }
    else if (strcmp(system, "figure-8") == 0) 
    {
        *objects_count = 3;
        *x = malloc(*objects_count * 3 * sizeof(real));
        *v = malloc(*objects_count * 3 * sizeof(real));
        *m = malloc(*objects_count * sizeof(real));

        (*x)[0] = 0.970043;
        (*x)[1] = -0.24308753;
        (*x)[2] = 0.0;

        (*x)[3] = -0.970043;
        (*x)[4] = 0.24308753;
        (*x)[5] = 0.0;

        (*x)[6] = 0.0;
        (*x)[7] = 0.0;
        (*x)[8] = 0.0;

        (*v)[0] = 0.466203685;
        (*v)[1] = 0.43236573;
        (*v)[2] = 0.0;

        (*v)[3] = 0.466203685;
        (*v)[4] = 0.43236573;
        (*v)[5] = 0.0;

        (*v)[6] = -0.93240737;
        (*v)[7] = -0.86473146;
        (*v)[8] = 0.0;

        (*m)[0] = 1.0 / *G;
        (*m)[1] = 1.0 / *G;
        (*m)[2] = 1.0 / *G;
    }
    else if (strcmp(system, "pyth-3-body") == 0) 
    {
        *objects_count = 3;
        *x = malloc(*objects_count * 3 * sizeof(real));
        *v = malloc(*objects_count * 3 * sizeof(real));
        *m = malloc(*objects_count * sizeof(real));

        (*x)[0] = 1.0;
        (*x)[1] = 3.0;
        (*x)[2] = 0.0;

        (*x)[3] = -2.0;
        (*x)[4] = -1.0;
        (*x)[5] = 0.0;

        (*x)[6] = 1.0;
        (*x)[7] = -1.0;
        (*x)[8] = 0.0;

        (*v)[0] = 0.0;
        (*v)[1] = 0.0;
        (*v)[2] = 0.0;

        (*v)[3] = 0.0;
        (*v)[4] = 0.0;
        (*v)[5] = 0.0;

        (*v)[6] = 0.0;
        (*v)[7] = 0.0;
        (*v)[8] = 0.0;

        (*m)[0] = 3.0 / *G;
        (*m)[1] = 4.0 / *G;
        (*m)[2] = 5.0 / *G;
    }
    else if (strcmp(system, "solar_system") == 0) 
    {
        *objects_count = 9;
        *x = malloc(*objects_count * 3 * sizeof(real));
        *v = malloc(*objects_count * 3 * sizeof(real));
        *m = malloc(*objects_count * sizeof(real));

        (*m)[0] = MASS_SUN;
        (*m)[1] = MASS_MERCURY;
        (*m)[2] = MASS_VENUS;
        (*m)[3] = MASS_EARTH;
        (*m)[4] = MASS_MARS;
        (*m)[5] = MASS_JUPITER;
        (*m)[6] = MASS_SATURN;
        (*m)[7] = MASS_URANUS;
        (*m)[8] = MASS_NEPTUNE;

        real R_CM[3];
        real V_CM[3];
        const real M = (
            (*m)[0] + (*m)[1] + (*m)[2] + (*m)[3] + (*m)[4] 
            + (*m)[5] + (*m)[6] + (*m)[7] + (*m)[8]
        );

        for (int i = 0; i < 3; i++)
        {
            R_CM[i] = 1 / M * (
                (*m)[0] * POS_SUN[i] + (*m)[1] * POS_MERCURY[i] + (*m)[2] * POS_VENUS[i]
                + (*m)[3] * POS_EARTH[i] + (*m)[4] * POS_MARS[i] + (*m)[5] * POS_JUPITER[i]
                + (*m)[6] * POS_SATURN[i] + (*m)[7] * POS_URANUS[i] + (*m)[8] * POS_NEPTUNE[i]
            );

            V_CM[i] = 1 / M * (
                (*m)[0] * VEL_SUN[i] + (*m)[1] * VEL_MERCURY[i] + (*m)[2] * VEL_VENUS[i]
                + (*m)[3] * VEL_EARTH[i] + (*m)[4] * VEL_MARS[i] + (*m)[5] * VEL_JUPITER[i]
                + (*m)[6] * VEL_SATURN[i] + (*m)[7] * VEL_URANUS[i] + (*m)[8] * VEL_NEPTUNE[i]
            );
        }

        for (int i = 0; i < 3; i++)
        {
            (*x)[0 * 3 + i] = POS_SUN[i] - R_CM[i];
            (*x)[1 * 3 + i] = POS_MERCURY[i] - R_CM[i];
            (*x)[2 * 3 + i] = POS_VENUS[i] - R_CM[i];
            (*x)[3 * 3 + i] = POS_EARTH[i] - R_CM[i];
            (*x)[4 * 3 + i] = POS_MARS[i] - R_CM[i];
            (*x)[5 * 3 + i] = POS_JUPITER[i] - R_CM[i];
            (*x)[6 * 3 + i] = POS_SATURN[i] - R_CM[i];
            (*x)[7 * 3 + i] = POS_URANUS[i] - R_CM[i];
            (*x)[8 * 3 + i] = POS_NEPTUNE[i] - R_CM[i];

            (*v)[0 * 3 + i] = VEL_SUN[i] - V_CM[i];
            (*v)[1 * 3 + i] = VEL_MERCURY[i] - V_CM[i];
            (*v)[2 * 3 + i] = VEL_VENUS[i] - V_CM[i];
            (*v)[3 * 3 + i] = VEL_EARTH[i] - V_CM[i];
            (*v)[4 * 3 + i] = VEL_MARS[i] - V_CM[i];
            (*v)[5 * 3 + i] = VEL_JUPITER[i] - V_CM[i];
            (*v)[6 * 3 + i] = VEL_SATURN[i] - V_CM[i];
            (*v)[7 * 3 + i] = VEL_URANUS[i] - V_CM[i];
            (*v)[8 * 3 + i] = VEL_NEPTUNE[i] - V_CM[i];
        }
    }
    else if (strcmp(system, "solar_system_plus") == 0) 
    {
        *objects_count = 12;
        *x = malloc(*objects_count * 3 * sizeof(real));
        *v = malloc(*objects_count * 3 * sizeof(real));
        *m = malloc(*objects_count * sizeof(real));
        
        (*m)[0] = MASS_SUN;
        (*m)[1] = MASS_MERCURY;
        (*m)[2] = MASS_VENUS;
        (*m)[3] = MASS_EARTH;
        (*m)[4] = MASS_MARS;
        (*m)[5] = MASS_JUPITER;
        (*m)[6] = MASS_SATURN;
        (*m)[7] = MASS_URANUS;
        (*m)[8] = MASS_NEPTUNE;
        (*m)[9] = MASS_PLUTO;
        (*m)[10] = MASS_CERES;
        (*m)[11] = MASS_VESTA;

        real R_CM[3];
        real V_CM[3];
        const real M = (
            (*m)[0] + (*m)[1] + (*m)[2] + (*m)[3] + (*m)[4] + (*m)[5] 
            + (*m)[6] + (*m)[7] + (*m)[8] + (*m)[9] + (*m)[10] + (*m)[11]
        );

        for (int i = 0; i < 3; i++)
        {
            R_CM[i] = 1 / M * (
                (*m)[0] * POS_SUN[i] + (*m)[1] * POS_MERCURY[i] + (*m)[2] * POS_VENUS[i]
                + (*m)[3] * POS_EARTH[i] + (*m)[4] * POS_MARS[i] + (*m)[5] * POS_JUPITER[i]
                + (*m)[6] * POS_SATURN[i] + (*m)[7] * POS_URANUS[i] + (*m)[8] * POS_NEPTUNE[i]
                + (*m)[9] * POS_PLUTO[i] + (*m)[10] * POS_CERES[i] + (*m)[11] * POS_VESTA[i]
            );

            V_CM[i] = 1 / M * (
                (*m)[0] * VEL_SUN[i] + (*m)[1] * VEL_MERCURY[i] + (*m)[2] * VEL_VENUS[i]
                + (*m)[3] * VEL_EARTH[i] + (*m)[4] * VEL_MARS[i] + (*m)[5] * VEL_JUPITER[i]
                + (*m)[6] * VEL_SATURN[i] + (*m)[7] * VEL_URANUS[i] + (*m)[8] * VEL_NEPTUNE[i]
                + (*m)[9] * VEL_PLUTO[i] + (*m)[10] * VEL_CERES[i] + (*m)[11] * VEL_VESTA [i]
            );
        }

        for (int i = 0; i < 3; i++)
        {
            (*x)[0 * 3 + i] = POS_SUN[i] - R_CM[i];
            (*x)[1 * 3 + i] = POS_MERCURY[i] - R_CM[i];
            (*x)[2 * 3 + i] = POS_VENUS[i] - R_CM[i];
            (*x)[3 * 3 + i] = POS_EARTH[i] - R_CM[i];
            (*x)[4 * 3 + i] = POS_MARS[i] - R_CM[i];
            (*x)[5 * 3 + i] = POS_JUPITER[i] - R_CM[i];
            (*x)[6 * 3 + i] = POS_SATURN[i] - R_CM[i];
            (*x)[7 * 3 + i] = POS_URANUS[i] - R_CM[i];
            (*x)[8 * 3 + i] = POS_NEPTUNE[i] - R_CM[i];
            (*x)[9 * 3 + i] = POS_PLUTO[i] - R_CM[i];
            (*x)[10 * 3 + i] = POS_CERES[i] - R_CM[i];
            (*x)[11 * 3 + i] = POS_VESTA[i] - R_CM[i];

            (*v)[0 * 3 + i] = VEL_SUN[i] - V_CM[i];
            (*v)[1 * 3 + i] = VEL_MERCURY[i] - V_CM[i];
            (*v)[2 * 3 + i] = VEL_VENUS[i] - V_CM[i];
            (*v)[3 * 3 + i] = VEL_EARTH[i] - V_CM[i];
            (*v)[4 * 3 + i] = VEL_MARS[i] - V_CM[i];
            (*v)[5 * 3 + i] = VEL_JUPITER[i] - V_CM[i];
            (*v)[6 * 3 + i] = VEL_SATURN[i] - V_CM[i];
            (*v)[7 * 3 + i] = VEL_URANUS[i] - V_CM[i];
            (*v)[8 * 3 + i] = VEL_NEPTUNE[i] - V_CM[i];
            (*v)[9 * 3 + i] = VEL_PLUTO[i] - V_CM[i];
            (*v)[10 * 3 + i] = VEL_CERES[i] - V_CM[i];
            (*v)[11 * 3 + i] = VEL_VESTA[i] - V_CM[i];
        }
    }
    else
    {
        *custom_sys_flag = 1;
    }
}

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
                for (int k = 0; k < 3; k++)
                {
                    temp_vec[k] = (
                        sol_state[*count][i * 3 + k]
                        - sol_state[*count][j * 3 + k]
                    );
                }
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

WIN32DLL_API void acceleration(int objects_count, real *restrict x, real *restrict a, const real *restrict m, real G)
{   
    real R_norm;
    real temp_value;
    real temp_vec[3];
    real R[3];

    // Empty the input array
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    for(int i = 0; i < objects_count; i++)
    {
        for(int j = i + 1; j < objects_count; j++)
        {
            // Calculate \vec{R} and its norm
            R[0] = x[i * 3 + 0] - x[j * 3 + 0];
            R[1] = x[i * 3 + 1] - x[j * 3 + 1];
            R[2] = x[i * 3 + 2] - x[j * 3 + 2];
            R_norm = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[i * 3 + 0] -= temp_vec[0] * m[j];
            a[i * 3 + 1] -= temp_vec[1] * m[j];
            a[i * 3 + 2] -= temp_vec[2] * m[j];
            a[j * 3 + 0] += temp_vec[0] * m[i];
            a[j * 3 + 1] += temp_vec[1] * m[i];
            a[j * 3 + 2] += temp_vec[2] * m[i];
        }
    }
}

WIN32DLL_API Solutions euler(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
)
{   
    // Initialize system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int custom_sys_flag = 0;
    initialize_system(system, &x, &v, &m, &objects_count, &G, &custom_sys_flag);

    // Custom system
    if (custom_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

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
    
    // Allocate memory for calculation
    real *temp_x = malloc(objects_count * 3 * sizeof(real));
    real *temp_v = malloc(objects_count * 3 * sizeof(real));
    real *a = malloc(objects_count * 3 * sizeof(real));
    
    // Allocate memory for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(store_npts * sizeof(double));
    double *sol_dt = malloc(store_npts * sizeof(double));

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

    // Main Loop
    while ((count + 1) <= npts)
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
        if ((count + 1) == npts)
        {
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sol_state[(store_npts - 1) * objects_count * 6 + i * 3 + j] = x[i * 3 + j];
                    sol_state[(store_npts - 1) * objects_count * 6 + objects_count * 3 + i * 3 + j] = v[i * 3 + j];
                }
            }
            sol_time[store_npts - 1] = dt * count;
        }
        else if (((count + 1) % store_every_n == 0))
        {
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sol_state[(*store_count + 1) * objects_count * 6 + i * 3 + j] = x[i * 3 + j];
                    sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3 + i * 3 + j] = v[i * 3 + j];
                }
            }
            sol_time[*store_count + 1] = dt * (count + 1);
            if (!((*store_count + 1) == (store_npts - 1)))
            {
                (*store_count)++;
            }
        }
        count++;
    }
    free(x);
    free(v);
    free(a);
    free(temp_x);
    free(temp_v);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    return (Solutions) {sol_state, sol_time, sol_dt, m, G, objects_count};
}

WIN32DLL_API Solutions euler_cromer(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
)
{   
    // Initialize system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int custom_sys_flag = 0;
    initialize_system(system, &x, &v, &m, &objects_count, &G, &custom_sys_flag);

    // Custom system
    if (custom_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

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

    // Allocate memory for calculation
    real *temp_x = malloc(objects_count * 3 * sizeof(real));
    real *temp_v = malloc(objects_count * 3 * sizeof(real));
    real *a = malloc(objects_count * 3 * sizeof(real));
    
    // Allocate memory for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(store_npts * sizeof(double));
    double *sol_dt = malloc(store_npts * sizeof(double));

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

    // Main Loop
    while ((count + 1) <= npts)
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
        if ((count + 1) == npts)
        {
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sol_state[(store_npts - 1) * objects_count * 6 + i * 3 + j] = x[i * 3 + j];
                    sol_state[(store_npts - 1) * objects_count * 6 + objects_count * 3 + i * 3 + j] = v[i * 3 + j];
                }
            }
            sol_time[store_npts - 1] = dt * count;
        }
        else if (((count + 1) % store_every_n == 0))
        {
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sol_state[(*store_count + 1) * objects_count * 6 + i * 3 + j] = x[i * 3 + j];
                    sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3 + i * 3 + j] = v[i * 3 + j];
                }
            }
            sol_time[*store_count + 1] = dt * (count + 1);
            if (!((*store_count + 1) == (store_npts - 1)))
            {
                (*store_count)++;
            }
        }
        count++;
    }
    free(x);
    free(v);
    free(a);
    free(temp_x);
    free(temp_v);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    return (Solutions) {sol_state, sol_time, sol_dt, m, G, objects_count};
}

WIN32DLL_API Solutions rk4(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
)
{
    // Initialize system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int custom_sys_flag = 0;
    initialize_system(system, &x, &v, &m, &objects_count, &G, &custom_sys_flag);

    // Custom system
    if (custom_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

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

    // Allocate memory for calculation
    real *temp_x = malloc(objects_count * 3 * sizeof(real));
    real *temp_v = malloc(objects_count * 3 * sizeof(real));
    real *a = malloc(objects_count * 3 * sizeof(real));
    real *vk1 = malloc(objects_count * 3 * sizeof(real));
    real *vk2 = malloc(objects_count * 3 * sizeof(real));
    real *vk3 = malloc(objects_count * 3 * sizeof(real));
    real *vk4 = malloc(objects_count * 3 * sizeof(real));
    real *xk1 = malloc(objects_count * 3 * sizeof(real));
    real *xk2 = malloc(objects_count * 3 * sizeof(real));
    real *xk3 = malloc(objects_count * 3 * sizeof(real));
    real *xk4 = malloc(objects_count * 3 * sizeof(real));

    // Allocate memory for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(store_npts * sizeof(double));
    double *sol_dt = malloc(store_npts * sizeof(double));

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

    // Main Loop
    while ((count + 1) <= npts)
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
        if ((count + 1) == npts)
        {
            for (int j = 0; j < objects_count; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    sol_state[(store_npts - 1) * objects_count * 6 + j * 3 + k] = x[j * 3 + k];
                    sol_state[(store_npts - 1) * objects_count * 6 + objects_count * 3 + j * 3 + k] = v[j * 3 + k];
                }
            }
            sol_time[store_npts - 1] = dt * count;
        }
        else if (((count + 1) % store_every_n == 0))
        {
            for (int j = 0; j < objects_count; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    sol_state[(*store_count + 1) * objects_count * 6 + j * 3 + k] = x[j * 3 + k];
                    sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3 + j * 3 + k] = v[j * 3 + k];
                }
            }
            sol_time[*store_count + 1] = dt * (count + 1);
            if (!((*store_count + 1) == (store_npts - 1)))
            {
                (*store_count)++;
            }
        }
        count++;
    }

    free(x);
    free(v);
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

    return (Solutions) {sol_state, sol_time, sol_dt, m, G, objects_count};
}

WIN32DLL_API Solutions leapfrog(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
)
{   
    // Initialize system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int custom_sys_flag = 0;
    initialize_system(system, &x, &v, &m, &objects_count, &G, &custom_sys_flag);

    // Custom system
    if (custom_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

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

    // Allocate memory for calculation
    real *temp_x = malloc(objects_count * 3 * sizeof(real));
    real *temp_v = malloc(objects_count * 3 * sizeof(real));
    real *a_0 = malloc(objects_count * 3 * sizeof(real));
    real *a_1 = malloc(objects_count * 3 * sizeof(real));

    acceleration(objects_count, x, a_1, m, G);

    // Allocate memory for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(store_npts * sizeof(double));
    double *sol_dt = malloc(store_npts * sizeof(double));

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

    // Main Loop
    while ((count + 1) <= npts)
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
        if ((count + 1) == npts)
        {
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sol_state[(store_npts - 1) * objects_count * 6 + i * 3 + j] = x[i * 3 + j];
                    sol_state[(store_npts - 1) * objects_count * 6 + objects_count * 3 + i * 3 + j] = v[i * 3 + j];
                }
            }
            sol_time[store_npts - 1] = dt * count;
        }
        else if (((count + 1) % store_every_n == 0))
        {
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    sol_state[(*store_count + 1) * objects_count * 6 + i * 3 + j] = x[i * 3 + j];
                    sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3 + i * 3 + j] = v[i * 3 + j];
                }
            }
            sol_time[*store_count + 1] = dt * (count + 1);
            if (!((*store_count + 1) == (store_npts - 1)))
            {
                (*store_count)++;
            }
        }
        count++;
    }

    free(x);
    free(v);
    free(temp_x);
    free(temp_v);
    free(a_0);
    free(a_1);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    return (Solutions) {sol_state, sol_time, sol_dt, m, G, objects_count};
}

WIN32DLL_API Solutions rk_embedded(
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
    int custom_sys_objects_count
)
{   
    // Initialize system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int custom_sys_flag = 0;
    initialize_system(system, &x, &v, &m, &objects_count, &G, &custom_sys_flag);

    // Custom system
    if (custom_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

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
    
    real *a = malloc(objects_count * 3 * sizeof(real));
    acceleration(objects_count, x, a, m, G);

    // Initialization
    int power;
    int power_test;
    real *coeff = NULL;
    int len_weights;
    real *weights = NULL;
    real *weights_test = NULL;
    rk_embedded_butcher_tableaus(order, &power, &power_test, &coeff, &len_weights, &weights, &weights_test);

    real abs_tolerance = input_abs_tolerance;
    real rel_tolerance = input_rel_tolerance;

    int stages = len_weights;
    int min_power = power < power_test ? power : power_test;

    real *error_estimation_delta_weights = malloc(len_weights * sizeof(real));
    for (int stage = 0; stage < stages; stage++)
    {
        error_estimation_delta_weights[stage] = weights[stage] - weights_test[stage];
    }

    // Safety factors for step-size control:
    real safety_fac_max = 6.0;
    real safety_fac_min = 0.33;
    real safety_fac = pow(0.38, (1.0 / (1.0 + (real) min_power)));

    // Initialize memory for calculation
    real sum, error, dt_new; 
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

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(NPTS * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(NPTS * sizeof(double));
    double *sol_dt = malloc(NPTS * sizeof(double));
    int buffer_size = NPTS;

    // Initial value
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sol_state[i * 3 + j] = x[i * 3 + j];
            sol_state[objects_count * 3 + i * 3 + j] = v[i * 3 + j];
        }
    }
    sol_time[*store_count + 1] = 0.0;
    real dt = rk_embedded_initial_dt(objects_count, power, x, v, a, m, G, abs_tolerance, rel_tolerance);
    sol_dt[*store_count + 1] = dt;

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
                for (int i = 0; i < objects_count; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        sol_state[(*store_count + 1) * objects_count * 6 + i * 3 + j] = x[i * 3 + j];
                        sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3 + i * 3 + j] = v[i * 3 + j];
                    }
                }
                *store_count += 1;
            }

            else if (*t >= tf)
            {
                sol_time[*store_count] = *t;
                sol_dt[*store_count] = dt;
                for (int i = 0; i < objects_count; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        sol_state[(*store_count + 1) * objects_count * 6 + i * 3 + j] = x[i * 3 + j];
                        sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3 + i * 3 + j] = v[i * 3 + j];
                    }
                } 
            }

            // End simulation as t >= tf
            if (*t >= tf)
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

                return (Solutions) {sol_state, sol_time, sol_dt, m, G, objects_count};
            }

            // Check buffer size and extend if full
            if ((*store_count + 1) == buffer_size)
            {   
                buffer_size += NPTS;
                sol_state = realloc(sol_state, buffer_size * objects_count * 6 * sizeof(real));
                sol_time = realloc(sol_time, buffer_size * sizeof(real));
                sol_dt = realloc(sol_dt, buffer_size * sizeof(real));
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
}

/* 
*   Calculate the initial time step for the Embedded RK method
*
*   Modified: Return dt * 1e-2 since this function gives initial dt thats too large
*/
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

/*
*   Butcher tableaus for embedded rk
*
*   Parameters: order (Must be one of 45 / 54 / 78 / 65)
*               power, power_test, coeff, weights, weights_test
*/
WIN32DLL_API void rk_embedded_butcher_tableaus(
    int order,
    int *restrict power,
    int *restrict power_test,
    real **coeff,
    int *restrict len_weights,
    real **weights,
    real **weights_test
)
{
    /*  Select integrator
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
            printf("Error: Unknown order for RK Embedded Integrator");
            exit(EXIT_FAILURE);
            break;
    }
}

WIN32DLL_API Solutions ias15(
    const char *restrict system,
    double *restrict t,
    double tf, 
    double input_tolerance,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count
)
{   
    // Initialize system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int custom_sys_flag = 0;
    initialize_system(system, &x, &v, &m, &objects_count, &G, &custom_sys_flag);

    // Custom system
    if (custom_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

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

    real *a = malloc(objects_count * 3 * sizeof(real));
    acceleration(objects_count, x, a, m, G);

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
    real *nodes = calloc(dim_nodes, sizeof(real));
    ias15_radau_spacing(nodes);
    real *aux_c = calloc(7 * 7, sizeof(real));
    ias15_aux_c(aux_c);
    real *aux_r = calloc(8 * 8, sizeof(real));
    ias15_aux_r(aux_r);
    real *aux_b0 = calloc((dim_nodes - 1) * objects_count * 3, sizeof(real));
    real *aux_b = calloc((dim_nodes - 1) * objects_count * 3, sizeof(real));
    real *aux_g = calloc((dim_nodes - 1) * objects_count * 3, sizeof(real));
    real *aux_e = calloc((dim_nodes - 1) * objects_count * 3, sizeof(real));

    int ias15_refine_flag = 0;

    // Arrays for ias15_step
    real *aux_a = calloc(dim_nodes * objects_count * 3, sizeof(real));
    real *x_step = calloc(objects_count * 3, sizeof(real));
    real *v_step = calloc(objects_count * 3, sizeof(real));
    real *a_step = calloc(objects_count * 3, sizeof(real));
    real *delta_b7 = calloc(objects_count * 3, sizeof(real));

    // Array for compute aux_g
    real *F = calloc(8 * objects_count * 3, sizeof(real));

    // Array for refine aux_b
    real *delta_aux_b = calloc(dim_nodes_minus_1 * objects_count * 3, sizeof(real));

    // Arrays for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *temp_x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *temp_v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(NPTS * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(NPTS * sizeof(double));
    double *sol_dt = malloc(NPTS * sizeof(double));
    int buffer_size = NPTS;

    // Initial value
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sol_state[i * 3 + j] = x[i * 3 + j];
            sol_state[objects_count * 3 + i * 3 + j] = v[i * 3 + j];
        }
    }
    sol_time[*store_count + 1] = 0.0;
    real dt = ias15_initial_dt(15, x, v, a, m, objects_count, G);
    sol_dt[*store_count + 1] = dt;

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
            &dt,
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
            temp_v_err_comp_sum
        );

        // Update count
        count += 1;

        // Store step
        if ((count + 1) % store_every_n == 0)
        {
            sol_time[*store_count + 1] = *t;
            sol_dt[*store_count + 1] = dt;
            for (int j = 0; j < objects_count; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    sol_state[(*store_count + 1) * objects_count * 6 + j * 3 + k] = x[j * 3 + k];
                    sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3 + j * 3 + k] = v[j * 3 + k];
                }
            }  
            *store_count += 1;
        }
        else if (*t >= tf)
        {
            sol_time[*store_count + 1] = *t;
            sol_dt[*store_count + 1] = dt;
            for (int j = 0; j < objects_count; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    sol_state[(*store_count + 1) * objects_count * 6 + j * 3 + k] = x[j * 3 + k];
                    sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3 + j * 3 + k] = v[j * 3 + k];
                }
            }  
        }

        // End simulation as t >= tf
        if (*t >= tf)
        {
            free(x);
            free(v);
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
            
            return (Solutions) {sol_state, sol_time, sol_dt, m, G, objects_count};
        }

        // Check buffer size and extend if full
        if ((*store_count + 1) == buffer_size)
        {   
            buffer_size += NPTS;
            sol_state = realloc(sol_state, buffer_size * objects_count * 6 * sizeof(real));
            sol_time = realloc(sol_time, buffer_size * sizeof(real));
            sol_dt = realloc(sol_dt, buffer_size * sizeof(real));
        }
    }
}

// Advance IAS15 for one step
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
    real *restrict temp_v_err_comp_sum
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
                acceleration(objects_count, x, &aux_a[i * objects_count * 3], m, G);

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
        acceleration(objects_count, x, a, m, G);

        // Estimate relative error
        error_b7 = abs_max_vec(&aux_b[dim_nodes_minus_2 * objects_count * 3], objects_count * 3) / abs_max_vec(a, objects_count * 3);
        error = pow((error_b7 / tolerance), exponent);

        // Step-size for the next step
        if (error != 0.0)
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

// Calculate position in the predictor-corrector algorithm to calculate aux_b and aux_g
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

// Calculate velocity in the predictor-corrector algorithm to calculate aux_b and aux_g
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

// Calculate the position of the next step
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

// Calculate the velocity of the next step
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

/*
*   Initialize the radau spacing nodes for IAS15
*
*   Return: pointer to 1D array with length 8 in data type real
*/
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

/*
*   Initialize the auxiliary coefficients aux_c for IAS15
*
*   Parameter: pointer to 1D array with length 7 * 7 in data type real
*/
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

/*
*   Initialize auxiliary coefficients aux_r for IAS15
*
*   Return: pointer to 1D array with length 8 * 8 in data type real
*/
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

/*
*   Calculate the initial time step for IAS15
*
*   return type: real
*/
real ias15_initial_dt(
    int power,
    real *restrict x,
    real *restrict v,
    real *restrict a,
    const real *m,
    int objects_count,
    real G
)
{
    real d_0 = abs_max_vec(x, objects_count * 3);
    real d_1 = abs_max_vec(a, objects_count * 3);
    real d_2;
    real dt_0;
    real dt_1;
    real *x_1 = malloc(objects_count * 3 * sizeof(real));
    real *a_1 = malloc(objects_count * 3 * sizeof(real));

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
    acceleration(objects_count, x_1, a_1, m, G);

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
