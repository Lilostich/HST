#include <iostream>

#include <vector>
#include <math.h>
#include <time.h>
#include <fstream>
#include <fcntl.h>

class Matrix {
public:
    Matrix(std::vector<double> matr);
    Matrix(std::string fileName);
    Matrix(){}
    void initSize(int n);

    double& getValue(int i, int j);

    std::vector<double> multiply(std::vector<double> rightVector,time_t *sum);

    void save(std::string fileName);
    int getSize(){return size;}
private:
    std::vector<double> values;
    int size;
};


Matrix::Matrix(std::string fileName) {
    FILE *fd;

    fd = fopen(fileName.c_str(),"r+b");
    fscanf(fd,"%d ",&size);
    double cell = 0;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            fscanf(fd,"%lf ",&cell);
            values.push_back(cell);
        }
    }

    fclose(fd);
}

///
/// \brief Умножение вектора на вектор
/// \param[in] matrix - матрица
/// \param[in] vect - вектор
/// \param[out] value - результат умножения
/// \param size - размер векторов
///
__global__ void kernelMultiply(double *matrix, double *vect, double *result,double* debug, int size){
    int tid = (threadIdx.x + blockIdx.x * blockDim.x)%size;
    int i = threadIdx.x;
    int j = blockIdx.x;
    result[i] += matrix[i * size + j] * vect[j];

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

Matrix::Matrix(std::vector<double> matr) :values(matr){size = sqrt(matr.size());}

void Matrix::initSize(int n) {values.reserve(n*n); size = n*n;}

double &Matrix::getValue(int i, int j) {return values[j + i * size];}

std::vector<double> Matrix::multiply(std::vector<double> rightVector,time_t *time) {
    std::vector<double> _values(values);
    std::vector<double> res;
//    printf("1\n");

    double *rVector = (double*)malloc(sizeof(double) * size);
    double *matrix = (double*)malloc(sizeof(double) * size * size);
    double *result = (double*)malloc(sizeof(double) * size);

    double *d_rVector, *d_matrix, *d_result;

//    for(int ii = 0; ii < res.size(); ii++){
//    }


    for(int i = 0; i < size; i++){
        rVector[i] = rightVector[i];
        result[i] = 0;
//        if(i < 10)
//            std::cout << rVector[i] << " ";
    }
//    std::cout << std::endl;
//    printf("\n");
    for(int i = 0; i < size * size; i++){
        matrix[i] = _values[i];
//        if(i < 10)
//            std::cout << matrix[i] << " ";
    }
//    std::cout << std::endl;

    cudaMalloc((void**)&d_rVector,sizeof(double) * size);
    cudaMalloc((void**)&d_matrix,sizeof(double) * size * size);
    cudaMalloc((void**)&d_result,sizeof(double) * size);
//    printf("asdasd2\n");

    cudaMemcpy(d_rVector,rVector,sizeof(double) * size,      cudaMemcpyHostToDevice);
    cudaMemcpy(d_matrix,matrix,sizeof(double) * size * size, cudaMemcpyHostToDevice);
//    printf("asdasd3\n");

//    cudaEvent_t start,stop;
    float gpuTime;
//    cudaEventCreate(&start);
//    cudaEventCreate(&stop);
//    cudaEventRecord(start, 0);
//    printf("asdasd4\n");
    double* debug = (double*)malloc(sizeof(double) * 10);
    double *d_debug;
    cudaMalloc((void**)&d_debug,sizeof(double) * 10);
    cudaMemcpy(d_debug,debug,sizeof(double) * 10, cudaMemcpyHostToDevice);
    time_t first = clock();
    kernelMultiply <<<size, size >>> (d_matrix, d_rVector, d_result, d_debug, size);
    time_t end = clock();

    cudaMemcpy(debug,d_debug,sizeof(double) * 10, cudaMemcpyDeviceToHost);
//    for(int j = 0; j < 10; j++){
//        std::cout << "deb " << j <<" "  << debug[j] << std::endl;
//    }
    cudaFree(d_debug);
    free(debug);

    *time += end-first;
//    printf("asdasd5\n");

//    cudaEventRecord(stop, 0);
//    cudaEventSynchronize(stop);
//    cudaEventElapsedTime(&gpuTime, start, stop);
//    printf("time on GPU = %.2lf ms \n", gpuTime);
//    cudaEventDestroy(start);
//    cudaEventDestroy(stop);
//    printf("asdasd6\n");

    cudaMemcpy(result, d_result, sizeof(double) * size, cudaMemcpyDeviceToHost);
//    printf("asdasd7\n");

    cudaFree(d_rVector);
    cudaFree(d_matrix);
    cudaFree(d_result);
//    printf("asdasd8\n");


    free(rVector);
    free(matrix);
//    printf("asdasd9\n");

    res.reserve(size);
    res.resize(size);
//    printf("asdasd10\n");
    for(int i = 0; i < size; i++){
        res[i] = result[i];
//        if(i < 10)
//            std::cout << result[i] << " " << res[i] << std::endl;
    }

//    for(int ii = 0; ii < res.size(); ii++){
//        if(ii < 10)
//            std::cout << result[ii] << " ";
//    }
//    std::cout << std::endl;

    free(result);
//    printf("asdasd11\n");
    return res;
}

void Matrix::save(std::string fileName) {
    std::ofstream fd;

    fd.open(fileName);
    fd << this->size << " ";
    for(int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            fd << values[i * size + j] << " ";
        }
    }
}

const double ACCURACY = 0.00001;
const std::string fileNames[] {"testMatrix.txt",
                               "matrixFile10.txt",
                               "matrixFile50.txt",
                               "matrixFile100.txt",
                               "matrixFile250.txt",
                               "matrixFile500.txt"};
const size_t FILE_COUNT = 6;
const size_t MAX_ITERATION = 100;

/// Длина вектора
double _dist(const std::vector<double> &rv){
    double res = 0;
    for (unsigned int i = 0; i < rv.size(); i++) {
        res += rv[i] * rv[i];
    }
    res = sqrt(res);
    return res != 0 ? res : 1;
}

/// деление вектора на значение
std::vector<double> _div(std::vector<double> rv, double diver){
    if (diver == 0)
        return rv;
    for (unsigned int i = 0; i < rv.size(); i++){
        rv[i] /= diver;
    }
    return rv;
}

void initRand(std::vector<double> &rv, int size){
//    printf("initRand\n");/

    rv.clear();
    rv.reserve(size);
    for (int i = 0; i < size; i++){
        rv.push_back(double(rand()));
    }
//    printf("endinitRand\n");

}

bool isVectorsNear(std::vector<double> &v1,std::vector<double> &v2){
//    printf("isVectorNear\n");
    for (int i = 0; i < v1.size(); i++){
        if (fabs(v1[i] - v2[i]) > ACCURACY)
        {
//            printf("endisVectorNear\n");
            return false;
        }
    }
//    printf("endisVectorNear\n");
    return true;
}

void findMyValue(Matrix &m,std::vector<double> &rv){
//    printf("findMyValue\n");
    std::vector<double> rvold;
    time_t sum = 0;
    int j = 0;
//    printf("%d\n",j++);
    time_t a = clock();
    time_t summ = 0;
    for (int i = 0; i < MAX_ITERATION; i++){
        rvold = rv;
//        printf("fmv1 %d\n",i);
//        for(int ii = 0; ii < rv.size(); ii++){
//            std::cout << rv[ii] << std::endl;
//        }
        std::vector<double> aa = m.multiply(rv,&summ);

        double distt = _dist(aa);
//        printf("fmv3 %d\n",i);
        rv = _div(aa, distt);
//        printf("fmv4 %d\n",i);
        if (isVectorsNear(rv,rvold))
            break;
//        printf("fmv5 %d\n",i);
    }
//    printf("%d\n",j++);

    time_t b = clock();
    sum += b - a;
    std::cout << "time is " << sum << std::endl;
    std::cout << "real time is " << summ << std::endl;
//    printf("endfindMyValue\n");
}

void saveVector(std::vector<double> &myVector, std::string fileName) {
//    printf("saveVector\n");
    FILE *fd;
    fd = fopen(("result_" + fileName).c_str(), "w+t");
    for (const auto &value : myVector){
        fprintf(fd, "%lf ", value);
    }
    fclose(fd);
//    printf("endSaveVector\n");
}

int main(int argc, char *argv[])
{
    for (int i = 0; i < FILE_COUNT; i++){
//        printf("iteration %d\n",i);
        Matrix matrix(fileNames[i]);
        int j = 0;
//        printf("%d\n",j++);
        std::vector<double> myVector;
        initRand(myVector,matrix.getSize());
//        printf("%d\n",j++);
        findMyValue(matrix,myVector);
//        printf("%d\n",j++);
        saveVector(myVector,fileNames[i]);
//        printf("%d\n",j++);
//        printf("enditeration %d\n",i);
    }
    return EXIT_SUCCESS;
}

/*
time is 166
time is 23
time is 111
time is 246
time is 688
time is 1274
 */