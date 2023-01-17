//
// Created by Dan on 23.12.2022.
//
#include <iostream>

#include <vector>
#include <math.h>
#include <time.h>
#include <fstream>
#include <fcntl.h>

///
/// \brief Умножение вектора на вектор
/// \param[in] matrix - матрица
/// \param[in] vect - вектор
/// \param[out] value - результат умножения
/// \param size - размер векторов
///
__global__ void kernelMultiply(double *lv, double *vect, double *result,double* debug, int size,int line){
    int i = threadIdx.x;
//    int i = blockIdx.x;
    result[line] += lv[i] * vect[i];

    debug[0] = threadIdx.x;
    debug[1] = blockIdx.x;
    debug[2] = blockDim.x;
    debug[3] = 5.6;
    debug[5] = tid;
    debug[6] = matrix[tid * size + 0];
    debug[7] = vect[0];
    debug[8] = size;
    debug[9] = matrix[6];

}

double multiply_str_str(std::vector<double> &leftVector, std::vector<double> &rightVector,time_t *time) {
    std::vector<double> _values(leftVector);
    std::vector<double> res;

    int size = rightVector.size();
    double *rVector = (double*)malloc(sizeof(double) * size);
    double *matrix = (double*)malloc(sizeof(double) * size);
    double *result = (double*)malloc(sizeof(double) * size);

    double *d_rVector, *d_matrix, *d_result;

    for(int i = 0; i < size; i++){
        rVector[i] = rightVector[i];
        result[i] = 0;
    }
    for(int i = 0; i < size * size; i++){
        matrix[i] = _values[i];
    }

    cudaMalloc((void**)&d_rVector,sizeof(double) * size);
    cudaMalloc((void**)&d_matrix,sizeof(double) * size * size);
    cudaMalloc((void**)&d_result,sizeof(double) * size);

    cudaMemcpy(d_rVector,rVector,sizeof(double) * size,      cudaMemcpyHostToDevice);
    cudaMemcpy(d_matrix,matrix,sizeof(double) * size * size, cudaMemcpyHostToDevice);

    double* debug = (double*)malloc(sizeof(double) * 10);
    double *d_debug;
    cudaMalloc((void**)&d_debug,sizeof(double) * 10);
    cudaMemcpy(d_debug,debug,sizeof(double) * 10, cudaMemcpyHostToDevice);
    time_t first = clock();
    kernelMultiply <<<size, size >>> (d_matrix, d_rVector, d_result, d_debug, size);
    time_t end = clock();

    cudaMemcpy(debug,d_debug,sizeof(double) * 10, cudaMemcpyDeviceToHost);
    cudaFree(d_debug);
    free(debug);

    *time += end-first;

    cudaMemcpy(result, d_result, sizeof(double) * size, cudaMemcpyDeviceToHost);

    cudaFree(d_rVector);
    cudaFree(d_matrix);
    cudaFree(d_result);

    free(rVector);
    free(matrix);

    res.reserve(size);
    res.resize(size);
    for(int i = 0; i < size; i++){
        res[i] = result[i];
    }

    free(result);
    return res;
}



