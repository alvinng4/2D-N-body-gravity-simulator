/**
 * \file integrator.h
 * \author Ching Yin Ng
 * \brief Function prototypes for all integrators
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "gravity_sim.h"

/**
 * \brief Euler first-order integrator
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param solutions Pointer to the solutions
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \retval SUCCESS If the simulation is successful
 * \retval error code If there is any error
 */
int euler(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
);

/**
 * \brief Euler-cromer integrator
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param solutions Pointer to the solutions
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \retval SUCCESS If the simulation is successful
 * \retval error code If there is any error
 */
int euler_cromer(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
);

/**
 * \brief Runge-Kutta 4th order integrator (RK4)
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param solutions Pointer to the solutions
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \retval SUCCESS If the simulation is successful
 * \retval error code If there is any error
 */
int rk4(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
);

/**
 * \brief LeapFrog integrator
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param solutions Pointer to the solutions
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \retval SUCCESS If the simulation is successful
 * \retval error code If there is any error
 */
int leapfrog(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
);

/**
 * \brief RK Embedded integrator
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param solutions Pointer to the solutions
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \retval SUCCESS If the simulation is successful
 * \retval error code If there is any error
 */
int rk_embedded(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
);

/**
 * \brief IAS15 integrator
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param solutions Pointer to the solutions
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \retval SUCCESS If the simulation is successful
 * \retval error code If there is any error
 */
int ias15(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
);

/**
 * \brief WHFast integrator
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param solutions Pointer to the solutions
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \retval SUCCESS If the simulation is successful
 * \retval error code If there is any error
 */
int whfast(
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