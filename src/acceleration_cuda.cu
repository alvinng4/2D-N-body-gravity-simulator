#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "acceleration_cuda.cuh"

WIN32DLL_API void acceleration_pairwise_cuda(
    int objects_count,
    double *__restrict x,
    double *__restrict a,
    const double *__restrict m,
    double G,
    double softening_length
)
{
    double *x_device = NULL;
    double *a_device = NULL;
    double *m_device = NULL;
    cudaError_t error = cudaSuccess;

    error = cudaMalloc((double **) &x_device, objects_count * 3 * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to allocate gpu memory for calculation\n");
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &a_device, objects_count * 3 * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to allocate gpu memory for calculation\n");
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &m_device, objects_count * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to allocate gpu memory for calculation\n");
        goto err_gpu_memory;
    }

    error = cudaMemcpy(x_device, x, objects_count * 3 * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to copy memory from host to device\n");
        goto err_gpu_memory;
    }
    error = cudaMemcpy(m_device, m, objects_count * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to copy memory from host to device\n");
        goto err_gpu_memory;
    }

    acceleration_pairwise_kernel <<< (objects_count + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>>(
        objects_count,
        x_device,
        a_device,
        m_device,
        G,
        softening_length
    );

    error = cudaMemcpy(a, a_device, objects_count * sizeof(double3), cudaMemcpyDeviceToHost);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to copy memory from device to host\n");
        goto err_gpu_memory;
    }

    cudaFree(x_device);
    cudaFree(a_device);
    cudaFree(m_device);
    return;

err_gpu_memory:
    cudaFree(x_device);
    cudaFree(a_device);
    cudaFree(m_device);
}

__global__ void acceleration_pairwise_kernel(
    int objects_count,
    double *__restrict x,
    double *__restrict a,
    const double *__restrict m,
    double G,
    double softening_length
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < objects_count)
    {
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;
        double obj_i_x = x[i * 3 + 0];
        double obj_i_y = x[i * 3 + 1];
        double obj_i_z = x[i * 3 + 2];

        for (int j = 0; j < objects_count; j++)
        {
            if (i != j)
            {
                double dx = x[j * 3 + 0] - obj_i_x;
                double dy = x[j * 3 + 1] - obj_i_y;
                double dz = x[j * 3 + 2] - obj_i_z;
                double r_norm = sqrt(dx * dx + dy * dy + dz * dz + softening_length * softening_length);

                double temp_value = G * m[j] / (r_norm * r_norm * r_norm);
                ax += temp_value * dx;
                ay += temp_value * dy;
                az += temp_value * dz;
            }
        }
        a[i * 3 + 0] = ax;
        a[i * 3 + 1] = ay;
        a[i * 3 + 2] = az;
    }
}

WIN32DLL_API void acceleration_pairwise_float_cuda(
    int objects_count,
    double *__restrict x,
    double *__restrict a,
    const double *__restrict m,
    double G,
    double softening_length
)
{
    float *x_device = NULL;
    float *a_device = NULL;
    float *m_device = NULL;
    float *x_float = NULL;
    float *a_float = NULL;
    float *m_float = NULL;
    cudaError_t error = cudaSuccess;

    error = cudaMalloc((float **) &x_device, objects_count * 3 * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to allocate gpu memory for calculation\n");
        goto err_gpu_memory;
    }
    error = cudaMalloc((float **) &a_device, objects_count * 3 * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to allocate gpu memory for calculation\n");
        goto err_gpu_memory;
    }
    error = cudaMalloc((float **) &m_device, objects_count * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to allocate gpu memory for calculation\n");
        goto err_gpu_memory;
    }

    x_float = (float *) malloc(objects_count * 3 * sizeof(float));
    a_float = (float *) malloc(objects_count * 3 * sizeof(float));
    m_float = (float *) malloc(objects_count * sizeof(float));
    if (x_float == NULL || a_float == NULL || m_float == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate cpu memory for calculation\n");
        goto err_cpu_memory;
    }

    for (int i = 0; i < objects_count; i++)
    {
        x_float[i * 3 + 0] = (float) x[i * 3 + 0];
        x_float[i * 3 + 1] = (float) x[i * 3 + 1];
        x_float[i * 3 + 2] = (float) x[i * 3 + 2];
        m_float[i] = (float) m[i];
    }

    error = cudaMemcpy(x_device, x_float, objects_count * 3 * sizeof(float), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to copy memory from host to device\n");
        goto err_gpu_memory;
    }
    error = cudaMemcpy(m_device, m_float, objects_count * sizeof(float), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to copy memory from host to device\n");
        goto err_gpu_memory;
    }

    acceleration_pairwise_float_kernel <<< (objects_count + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>>(
        objects_count,
        x_device,
        a_device,
        m_device,
        G,
        softening_length
    );

    error = cudaMemcpy(a_float, a_device, objects_count * 3 * sizeof(float), cudaMemcpyDeviceToHost);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to copy memory from device to host\n");
        goto err_gpu_memory;
    }
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = (double) a_float[i * 3 + 0];
        a[i * 3 + 1] = (double) a_float[i * 3 + 1];
        a[i * 3 + 2] = (double) a_float[i * 3 + 2];
    }

    free(x_float);
    free(a_float);
    free(m_float);
    cudaFree(x_device);
    cudaFree(a_device);
    cudaFree(m_device);
    return;

err_cpu_memory:
    free(x_float);
    free(a_float);
    free(m_float);
err_gpu_memory:
    cudaFree(x_device);
    cudaFree(a_device);
    cudaFree(m_device);
}

__global__ void acceleration_pairwise_float_kernel(
    int objects_count,
    float *__restrict x,
    float *__restrict a,
    const float *__restrict m,
    float G,
    float softening_length
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < objects_count)
    {
        float ax = 0.0;
        float ay = 0.0;
        float az = 0.0;
        float obj_i_x = x[i * 3 + 0];
        float obj_i_y = x[i * 3 + 1];
        float obj_i_z = x[i * 3 + 2];

        for (int j = 0; j < objects_count; j++)
        {
            if (i != j)
            {
                float dx = x[j * 3 + 0] - obj_i_x;
                float dy = x[j * 3 + 1] - obj_i_y;
                float dz = x[j * 3 + 2] - obj_i_z;
                float r_norm = sqrtf(dx * dx + dy * dy + dz * dz + softening_length * softening_length);

                float temp_value = G * m[j] / (r_norm * r_norm * r_norm);
                ax += temp_value * dx;
                ay += temp_value * dy;
                az += temp_value * dz;
            }
        }
        a[i * 3 + 0] = ax;
        a[i * 3 + 1] = ay;
        a[i * 3 + 2] = az;
    }
}

WIN32DLL_API void acceleration_pairwise_float_comp_sum_cuda(
    int objects_count,
    double *__restrict x,
    double *__restrict a,
    const double *__restrict m,
    double G,
    double softening_length
)
{
    float *x_device = NULL;
    float *a_device = NULL;
    float *m_device = NULL;
    float *x_float = NULL;
    float *a_float = NULL;
    float *m_float = NULL;
    cudaError_t error = cudaSuccess;

    error = cudaMalloc((float **) &x_device, objects_count * 3 * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to allocate gpu memory for calculation\n");
        goto err_gpu_memory;
    }
    error = cudaMalloc((float **) &a_device, objects_count * 3 * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to allocate gpu memory for calculation\n");
        goto err_gpu_memory;
    }
    error = cudaMalloc((float **) &m_device, objects_count * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to allocate gpu memory for calculation\n");
        goto err_gpu_memory;
    }

    x_float = (float *) malloc(objects_count * 3 * sizeof(float));
    a_float = (float *) malloc(objects_count * 3 * sizeof(float));
    m_float = (float *) malloc(objects_count * sizeof(float));
    if (x_float == NULL || a_float == NULL || m_float == NULL)
    {
        fprintf(stderr, "Error: Failed to allocate cpu memory for calculation\n");
        goto err_cpu_memory;
    }

    for (int i = 0; i < objects_count; i++)
    {
        x_float[i * 3 + 0] = (float) x[i * 3 + 0];
        x_float[i * 3 + 1] = (float) x[i * 3 + 1];
        x_float[i * 3 + 2] = (float) x[i * 3 + 2];
        m_float[i] = (float) m[i];
    }

    error = cudaMemcpy(x_device, x_float, objects_count * 3 * sizeof(float), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to copy memory from host to device\n");
        goto err_gpu_memory;
    }
    error = cudaMemcpy(m_device, m_float, objects_count * sizeof(float), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to copy memory from host to device\n");
        goto err_gpu_memory;
    }

    acceleration_pairwise_float_comp_sum_kernel <<< (objects_count + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>>(
        objects_count,
        x_device,
        a_device,
        m_device,
        G,
        softening_length
    );

    error = cudaMemcpy(a_float, a_device, objects_count * 3 * sizeof(float), cudaMemcpyDeviceToHost);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Error: Failed to copy memory from device to host\n");
        goto err_gpu_memory;
    }
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = (double) a_float[i * 3 + 0];
        a[i * 3 + 1] = (double) a_float[i * 3 + 1];
        a[i * 3 + 2] = (double) a_float[i * 3 + 2];
    }

    free(x_float);
    free(a_float);
    free(m_float);
    cudaFree(x_device);
    cudaFree(a_device);
    cudaFree(m_device);
    return;

err_cpu_memory:
    free(x_float);
    free(a_float);
    free(m_float);
err_gpu_memory:
    cudaFree(x_device);
    cudaFree(a_device);
    cudaFree(m_device);
}

__global__ void acceleration_pairwise_float_comp_sum_kernel(
    int objects_count,
    float *__restrict x,
    float *__restrict a,
    const float *__restrict m,
    float G,
    float softening_length
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < objects_count)
    {
        float ax = 0.0;
        float ay = 0.0;
        float az = 0.0;
        float obj_i_x = x[i * 3 + 0];
        float obj_i_y = x[i * 3 + 1];
        float obj_i_z = x[i * 3 + 2];
        float ax0 = 0.0;
        float ay0 = 0.0;
        float az0 = 0.0;
        float x_err_comp_sum = 0.0;
        float y_err_comp_sum = 0.0;
        float z_err_comp_sum = 0.0;

        for (int j = 0; j < objects_count; j++)
        {
            if (i != j)
            {
                float dx = x[j * 3 + 0] - obj_i_x;
                float dy = x[j * 3 + 1] - obj_i_y;
                float dz = x[j * 3 + 2] - obj_i_z;
                float r_norm = sqrtf(dx * dx + dy * dy + dz * dz + softening_length * softening_length);

                float temp_value = G * m[j] / (r_norm * r_norm * r_norm);

                ax0 = ax;
                ay0 = ay;
                az0 = az;

                x_err_comp_sum += temp_value * dx;
                y_err_comp_sum += temp_value * dy;
                z_err_comp_sum += temp_value * dz;

                ax = ax0 + x_err_comp_sum;
                ay = ay0 + y_err_comp_sum;
                az = az0 + z_err_comp_sum;

                x_err_comp_sum += ax0 - ax;
                y_err_comp_sum += ay0 - ay;
                z_err_comp_sum += az0 - az;
            }
        }
        a[i * 3 + 0] = ax;
        a[i * 3 + 1] = ay;
        a[i * 3 + 2] = az;
    }
}