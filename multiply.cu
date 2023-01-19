//
// Created by Dan on 23.12.2022.
//
// #include <iostream>

// #include <vector>
// #include <math.h>
// #include <time.h>
#include <stdio.h>
// #include <fcntl.h>
// #include <stdlib.h>

///
/// \brief Умножение вектора на вектор
/// \param[in] lv - левый вектор
/// \param[in] vect - правый вектор
/// \param[in] size - размер векторов
/// \param[out] result - результат умножения
/// \param[out] debug - вектор для дебаг значений


// __global__ void kernelMultiply(double *lv, double *vect, double *result, int size, double* debug){
__global__ void kernelMultiply(double *lv, double *vect, double *result,int size, double* debug){
    // result[0] = 404.0;
    // lv[0] = 404.0;
    // vect[0] = 404.0;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size)
        *result += lv[i] * vect[i];
    // if (debug[0] < (double)(int)threadIdx.x)
        // debug[0] = threadIdx.x;
    // if (debug[1] < (double)(int)blockIdx.x)
        // debug[1] = blockIdx.x;
    // if (debug[2] < (double)(int)blockDim.x)
        // debug[2] = blockDim.x;
}



///
/// \brief Умножение вектора на вектор
/// \param[in] matrix - матрица
/// \param[in] vect - вектор
/// \param[out] value - результат умножения
/// \param size - размер векторов
///
double multiply_str_str(double *leftVector, double *rightVector,int size, time_t *time) {
    
    cudaError_t error;
    int step = 0;
    double *d_rVector = NULL, *d_lVector = NULL, *d_result = NULL;
    double *result = (double *)malloc(sizeof(double));
    result[0] = -1;

    error = cudaMalloc((void**)&d_rVector,sizeof(double) * size);
    // printf("step %d error is %d\n",step++,error);
    error = cudaMalloc((void**)&d_lVector,sizeof(double) * size);
    // printf("step %d error is %d\n",step++,error);
    error = cudaMalloc((void**)&d_result,sizeof(double) * 1);
    // printf("step %d error is %d\n",step++,error);

    // printf("leftVector is [%lf, %lf, ...]\n",leftVector[0],leftVector[1]);
    // printf("rightVector is [%lf, %lf, ...]\n",rightVector[0],rightVector[1]);

    error = cudaMemcpy(d_rVector,rightVector,sizeof(double) * size, cudaMemcpyHostToDevice);
    // printf("step %d error is %d\n",step++,error);
    error = cudaMemcpy(d_lVector,leftVector ,sizeof(double) * size, cudaMemcpyHostToDevice);
    // printf("step %d error is %d\n",step++,error);
    error = cudaMemcpy(d_result,result ,sizeof(double), cudaMemcpyHostToDevice);
    // printf("step %d error is %d\n",step++,error);

    // DEBUG
    int debug_size = 10;
    double* debug = (double*)malloc(sizeof(double) * debug_size);
    debug[0] = debug[1] = debug[2] = -1;
    double *d_debug;
    error = cudaMalloc((void**)&d_debug,sizeof(double) * debug_size);
    // printf("step %d error is %d\n",step++,error);
    // \DEBUG

    error = cudaMemcpy(d_debug,debug,sizeof(double) * debug_size, cudaMemcpyHostToDevice);
    // printf("step %d error is %d\n",step++,error);
    time_t first = clock();
 
    // TRY RUN
    // kernelMultiply <<< 256, 256 >>> (d_lVector, d_rVector, d_result, size, d_debug);
    kernelMultiply <<< 1,1 >>> (d_lVector, d_rVector, d_result, size, d_debug);
    // printf("try run is is %d\n",error);
    cudaDeviceSynchronize();
    cudaThreadSynchronize(); 
    time_t end = clock();

    error = cudaMemcpy(debug,d_debug,sizeof(double) * debug_size, cudaMemcpyDeviceToHost);
    // printf("step %d error is %d\n",step++,error);
    // printf("max threadIdx.x is %lf\n", debug[0]);
    // printf("max blockIdx.x is %lf\n", debug[1]);
    // printf("max blockDim.x is %lf\n", debug[2]);
    error = cudaFree(d_debug);
    // printf("step %d error is %d\n",step++,error);
    // Вывод необходимых debug переменных
    free(debug);

    *time += end-first;
    
    error = cudaMemcpyAsync(result, d_result, sizeof(double), cudaMemcpyDeviceToHost);
    // printf("step %d error iss %d\n",step++,error);
    // printf("result is %lf\n",result[0]);
    

    error = cudaFree(d_rVector);
    // printf("step %d error is %d\n",step++,error);
    error = cudaFree(d_lVector);
    // printf("step %d error is %d\n",step++,error);
    error = cudaFree(d_result);
    // printf("step %d error is %d\n",step++,error);
    // printf("All steps ends\n");

    return *result;
}



