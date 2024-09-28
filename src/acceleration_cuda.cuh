#ifndef ACCELERATION_CUDA_H
#define ACCELERATION_CUDA_H

/**
 * acceleration_cuda: Acceleration calculation with CUDA
 */

#include <cuda_runtime.h>

#ifdef WIN32DLL_EXPORTS
    #define WIN32DLL_API __declspec(dllexport)
#else
    #define WIN32DLL_API 
#endif

#define BLOCK_SIZE 64

#ifdef __cplusplus
    extern "C" {
#endif
    WIN32DLL_API void acceleration_pairwise_cuda(
        int objects_count,
        double *__restrict x,
        double *__restrict a,
        const double *__restrict m,
        double G,
        double softening_length
    );

    WIN32DLL_API void acceleration_pairwise_float_cuda(
        int objects_count,
        double *__restrict x,
        double *__restrict a,
        const double *__restrict m,
        double G,
        double softening_length
    );

    WIN32DLL_API void acceleration_pairwise_float_comp_sum_cuda(
        int objects_count,
        double *__restrict x,
        double *__restrict a,
        const double *__restrict m,
        double G,
        double softening_length
    );
#ifdef __cplusplus
    }
#endif

__global__ void acceleration_pairwise_kernel(
    int objects_count,
    double *__restrict x,
    double *__restrict a,
    const double *__restrict m,
    double G,
    double softening_length
);

__global__ void acceleration_pairwise_float_kernel(
    int objects_count,
    float *__restrict x,
    float *__restrict a,
    const float *__restrict m,
    float G,
    float softening_length
);

__global__ void acceleration_pairwise_float_comp_sum_kernel(
    int objects_count,
    float *__restrict x,
    float *__restrict a,
    const float *__restrict m,
    float G,
    float softening_length
);

#endif