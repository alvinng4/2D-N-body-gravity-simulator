/**
 * \file storing.h
 * \author Ching Yin Ng
 * \brief Contains methods for storing the simulation output
 */

#ifndef STORING_H
#define STORING_H

#define STORING_METHOD_DEFAULT 1
#define STORING_METHOD_FLUSH 2
#define STORING_METHOD_DISABLED 3

/**
 * \brief Return storing method flag based on the input string
 * 
 * \param storing_method Name of the storing method
 * 
 * \retval SUCCESS If the storing method is recognized
 * \retval ERROR_UNKNOWN_STORING_METHOD If the storing method is not recognized
 */
int get_storing_method_flag(
    const char *restrict storing_method,
    uint *restrict storing_method_flag
);

/**
 * \brief Allocate memory for solutions output
 * 
 * \param solutions Pointer to the solutions
 * \param storing_param Pointer to the storing parameters
 * \param objects_count Number of objects in the system
 */
int allocate_solutions_memory(
    Solutions *restrict solutions,
    StoringParam *restrict storing_param,
    int objects_count
);

/**
 * \brief Open file for flushing
 * 
 * \param storing_param Pointer to the storing parameters
 * 
 * \retval SUCCESS If the file is opened successfully
 * \retval ERROR_FLUSH_FILE_OPEN If the file failed to open
 */
int open_flush_file(
    StoringParam *restrict storing_param
);

/**
 * \brief Close file for flushing
 * 
 * \param storing_param Pointer to the storing parameters
 * 
 * \retval SUCCESS If the file is closed successfully
 * \retval ERROR_FLUSH_FILE_CLOSE_NOT_FLUSH_METHOD If the storing method is not flush
 * \retval ERROR_FLUSH_FILE_CLOSE_IS_NULL If the file pointer is NULL
 * \retval ERROR_FLUSH_FILE_CLOSE If the file failed to close
 */
int close_flush_file(
    StoringParam *restrict storing_param
);

/**
 * \brief Flush solution step to CSV file
 * 
 * \param file Pointer to the file to store the solution
 * \param system Pointer to the gravitational system
 * \param simulation_status Pointer to the simulation status
 * \param solutions Pointer to the solutions
 * 
 * \retval SUCCESS If the solution is flushed successfully
 * \retval error code If there is any error
 */
int flush_solution_step_to_csv_file(
    FILE *restrict file,
    const System *restrict system,
    const SimulationStatus *restrict simulation_status,
    Solutions *restrict solutions
);

/**
 * \brief Store solution step to memory
 * 
 * \param system Pointer to the gravitational system
 * \param simulation_status Pointer to the simulation status
 * \param solutions Pointer to the solutions
 * 
 * \retval SUCCESS If the solution is stored successfully
 */
int store_solution_step_to_memory(
    const System *restrict system,
    const SimulationStatus *restrict simulation_status,
    Solutions *restrict solutions
);

/**
 * \brief Store solution step based on the storing parameters
 * 
 * \param storing_param Pointer to the storing parameters
 * \param system Pointer to the gravitational system
 * \param simulation_status Pointer to the simulation status
 * \param solutions Pointer to the solutions
 * 
 * \retval SUCCESS If the solution is stored successfully
 * \retval error code If there is any error
 */
int store_solution_step(
    StoringParam *restrict storing_param,
    const System *restrict system,
    const SimulationStatus *restrict simulation_status,
    Solutions *restrict solutions
);

/**
 * \brief Extend memory buffer for solution output
 * 
 * \param solutions Pointer to the solutions
 * \param storing_param Pointer to the storing parameters
 * \param objects_count Number of objects in the system
 * 
 * \retval SUCCESS If the memory buffer is extended successfully
 * \retval ERROR_SOL_OUTPUT_EXTEND_MEMORY_REALLOC If the memory reallocation failed
 */
int extend_sol_memory_buffer(
    Solutions *restrict solutions,
    StoringParam *restrict storing_param,
    const int objects_count
);

#endif