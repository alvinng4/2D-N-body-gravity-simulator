#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

WIN32DLL_API void euler(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    const real *restrict m, 
    real G, 
    real dt,
    int time_speed
)
{   
    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    // Main Loop
    for(int count = 0; count < time_speed; count++)
    {   
        acceleration(objects_count, x, a, m, G);
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            x[j][0] += v[j][0] * dt;
            x[j][1] += v[j][1] * dt;
            x[j][2] += v[j][2] * dt;
            v[j][0] += a[j][0] * dt;
            v[j][1] += a[j][1] * dt;
            v[j][2] += a[j][2] * dt;
        }    
    }

    free(a);
}

WIN32DLL_API void euler_cromer(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    const real *restrict m, 
    real G, 
    real dt,
    int time_speed
)
{   
    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    // Main Loop
    for(int count = 0; count < time_speed; count++)
    {   
        acceleration(objects_count, x, a, m, G);
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            v[j][0] += a[j][0] * dt;
            v[j][1] += a[j][1] * dt;
            v[j][2] += a[j][2] * dt;
            x[j][0] += v[j][0] * dt;
            x[j][1] += v[j][1] * dt;
            x[j][2] += v[j][2] * dt;
        }    
    }
    
    free(a);
}

WIN32DLL_API void rk4(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    const real *restrict m, 
    real G, 
    real dt,
    int time_speed
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

    // Main Loop
    for(int count = 0; count < time_speed; count++)
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


        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            v[j][0] += (vk1[j][0] + 2 * vk2[j][0] + 2 * vk3[j][0] + vk4[j][0]) * dt / 6.0;
            v[j][1] += (vk1[j][1] + 2 * vk2[j][1] + 2 * vk3[j][1] + vk4[j][1]) * dt / 6.0;
            v[j][2] += (vk1[j][2] + 2 * vk2[j][2] + 2 * vk3[j][2] + vk4[j][2]) * dt / 6.0;
            x[j][0] += (xk1[j][0] + 2 * xk2[j][0] + 2 * xk3[j][0] + xk4[j][0]) * dt / 6.0;
            x[j][1] += (xk1[j][1] + 2 * xk2[j][1] + 2 * xk3[j][1] + xk4[j][1]) * dt / 6.0;
            x[j][2] += (xk1[j][2] + 2 * xk2[j][2] + 2 * xk3[j][2] + xk4[j][2]) * dt / 6.0;
        }    
    } 

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

}

WIN32DLL_API void leapfrog(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real (*restrict a)[3], 
    const real *restrict m, 
    real G, 
    real dt,
    int time_speed
)
{   
    real (*a_0)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a_1)[3] = malloc(objects_count * 3 * sizeof(real));

    memcpy(a_1, a, objects_count * 3 * sizeof(real));

    // Main Loop
    for(int count = 0; count < time_speed; count++)
    {       
        // Use a_1 from last iteration as a_0
        memcpy(a_0, a_1, objects_count * 3 * sizeof(real));
        
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            x[j][0] += v[j][0] * dt + 0.5 * a_0[j][0] * dt * dt;
            x[j][1] += v[j][1] * dt + 0.5 * a_0[j][1] * dt * dt;
            x[j][2] += v[j][2] * dt + 0.5 * a_0[j][2] * dt * dt;
        }    
        acceleration(objects_count, x, a_1, m, G);
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            v[j][0] += 0.5 * (a_0[j][0] + a_1[j][0]) * dt;
            v[j][1] += 0.5 * (a_0[j][1] + a_1[j][1]) * dt;
            v[j][2] += 0.5 * (a_0[j][2] + a_1[j][2]) * dt;
        }    
    }

    memcpy(a, a_1, objects_count * 3 * sizeof(real));

    free(a_0);
    free(a_1);
}