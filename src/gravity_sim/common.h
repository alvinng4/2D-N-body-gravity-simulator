#ifndef COMMON_H
#define COMMON_H

/**
 * common: Contains commonly used functions for gravity 
 *         simulation, e.g. acceleration, vec_norm, etc.
 *         Definitions are also included in this file.
 */

#ifdef WIN32DLL_EXPORTS
    #define WIN32DLL_API __declspec(dllexport)
#else
    #define WIN32DLL_API 
#endif

typedef double real;


real abs_max_vec(const real *restrict vec, int vec_length);

real abs_max_vec_array(real (*restrict arr)[3], int objects_count);

real vec_norm(const real *restrict vec, int vec_length);

real compute_energy(
    int objects_count, 
    const real (*restrict x)[3],
    const real (*restrict v)[3],
    const real *restrict m, 
    real G
);

void acceleration(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict a)[3], 
    const real *restrict m, 
    real G
);


#endif