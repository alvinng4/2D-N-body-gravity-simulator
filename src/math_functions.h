#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include "gravity_sim.h"

/**
 * \brief Calculate the norm of a 3D vector
 * 
 * \param vec 1D array of size 3
 * 
 * \return Norm of the vector
 */
real vec_norm_3d(const real *restrict vec);

/**
 * \brief Calculate the norm of a vector
 * 
 * \param vec 1D array of size vec_length
 * \param vec_length Length of the vector
 * 
 * \return Norm of the vector
 */
real vec_norm(const real *restrict vec, int vec_length);

/**
 * \brief Calculate the absolute maximum value from a vector
 * 
 * \param vec 1D array of size vec_length
 * \param vec_length Length of the vector
 * 
 * \return Absolute maximum value from the vector
 */
real abs_max_vec(const real *restrict vec, int vec_length);

/**
 * \brief Calculate the dot product of two 3D vectors
 * 
 * \param vec_1 1D array of size 3
 * \param vec_2 1D array of size 3
 * 
 * \return Dot product of the two vectors
 */
real vec_dot_3d(
    const real *restrict vec_1,
    const real *restrict vec_2
);

/**
 * \brief Calculate the dot product of two vectors
 * 
 * \param vec_1 1D array of size vec_length
 * \param vec_2 1D array of size vec_length
 * \param vec_length Length of the vectors
 * 
 * \return Dot product of the two vectors
 */
real vec_dot(
    const real *restrict vec_1,
    const real *restrict vec_2,
    const int vec_length
);

/**
 * \brief Calculate 2^n using bit shifting
 * 
 * \param n Exponent
 * 
 * \return 2^n
 */
int fast_pow_of_2(int n);

#endif
