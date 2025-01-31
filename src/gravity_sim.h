/**
 * \file gravity_sim.h
 * \author Ching Yin Ng
 * \brief Contains important functions, structs, and definitions 
 *        for the gravity_sim library
 */

#ifndef GRAVITY_SIM_H
#define GRAVITY_SIM_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

// For exporting functions in Windows DLL as a dynamic-link library
#ifdef WIN32DLL_EXPORTS
    #define WIN32DLL_API __declspec(dllexport)
#else
    #define WIN32DLL_API 
#endif

#define IN_FILE static
#define ADAPTIVE_STEP_SIZE_SOL_BUFFER_SIZE 50000

typedef unsigned int uint;
typedef int64_t int64;
typedef double real;

typedef struct System
{
    real *x;
    real *v;
    real *m;
    int objects_count;
    real G;
} System;

typedef struct IntegratorParam
{
    const char *integrator;
    real dt;
    real tolerance;
    real initial_dt;
    real whfast_kepler_tol;
    int whfast_kepler_max_iter;
    bool whfast_kepler_auto_remove;
    real whfast_kepler_auto_remove_tol;
} IntegratorParam;

typedef struct AccelerationParam
{
    const char *method;
    real opening_angle;
    real softening_length;
    int order;
    uint acceleration_method_flag_;
} AccelerationParam;

typedef struct StoringParam
{
    const char *method;
    const char *flush_path;
    int storing_freq;
    uint storing_method_flag_;
    FILE *flush_file_;
    int64 max_sol_size_;
} StoringParam;

typedef struct Solutions 
{
    double *restrict sol_state;
    double *restrict sol_time;
    double *restrict sol_dt;
    int64 *sol_size_;
} Solutions;

// Current status of the simulation
typedef struct SimulationStatus
{
    real *t;
    real dt;
    real run_time_;
} SimulationStatus;

typedef struct Settings
{
    int verbose;
    bool *restrict is_exit;
} Settings;

typedef struct SimulationParam
{
    real tf;
    int64 n_steps_;
} SimulationParam;

/**
 * \brief Launch simulation from Python
 * 
 * \param x Pointer to the position array
 * \param v Pointer to the velocity array
 * \param m Pointer to the mass array
 * \param objects_count Number of objects in the system
 * \param G Gravitational constant
 * \param integrator Name of the integrator
 * \param dt Time step size
 * \param tolerance Tolerance for adaptive step size integrators
 * \param acceleration_method Name of the acceleration method
 * \param opening_angle Opening angle for the acceleration calculation
 * \param softening_length Softening length for the force calculation
 * \param order Order of the acceleration approximation
 * \param storing_method Name of the storing method
 * \param flush_path Path to the file to store the solution
 * \param storing_freq Storing frequency
 * \param sol_state Pointer of pointer to the solution state array to be updated
 * \param sol_time Pointer of pointer to the solution time array to be updated
 * \param sol_dt Pointer of pointer to the solution step size array to be updated
 * \param sol_size_ Pointer to the solution size to be updated
 * \param t Pointer to the current simulation time to be updated
 * \param simulation_status_last_dt Pointer to the last time step size of the simulation
 * \param run_time_ Pointer to the run time of the simulation to be updated
 * \param verbose Verbosity level
 * \param is_exit Pointer to the exit flag
 * \param tf Simulation time
 */
int launch_simulation_python(
    real *x,
    real *v,
    real *m,
    int objects_count,
    real G,
    const char *integrator,
    real dt,
    real tolerance,
    real initial_dt,
    real whfast_kepler_tol,
    int whfast_kepler_max_iter,
    bool whfast_kepler_auto_remove,
    real whfast_kepler_auto_remove_tol,
    const char *acceleration_method,
    real opening_angle,
    real softening_length,
    int order,
    const char *storing_method,
    const char *flush_path,
    int storing_freq,
    double **sol_state,
    double **sol_time,
    double **sol_dt,
    int64 *sol_size_,
    real *t,
    real *simulation_status_last_dt,
    real *run_time_,
    int verbose,
    bool *is_exit,
    real tf
);

/**
 * \brief Launch simulation
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param solutions Pointer to the solutions
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param tf Simulation time
 * 
 * \retval SUCCESS If the simulation is successful
 * \retval ERROR_SIMULATION_FAILURE If the simulation fails
 * 
 * \note Side effects: prints error message to stderr, unless settings->verbose <= 0
 */
int launch_simulation(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
);

#endif
