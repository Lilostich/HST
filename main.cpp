#include <iostream>

#include "Matrix.h"

const double ACCURACY = 0.00001;
const std::string fileNames[] {"matrixFile10.txt",
                         "matrixFile50.txt",
                         "matrixFile100.txt",
                         "matrixFile250.txt",
                         "matrixFile500.txt"};
const size_t FILE_COUNT = 5;
const size_t MAX_ITERATION = 100;

const int sizeMatrix[] {1144,2560,3620, 5724, 8095};

#define CHOICE 4

/// Длина вектора
double _dist(const std::vector<double> &rv){
    double res = 0;
    for (unsigned int i = 0; i < rv.size(); i++) {
        res += rv[i] * rv[i];
    }
    return sqrt(res);
}

/// деление вектора на значение
std::vector<double> _div(std::vector<double> rv, double diver){
    for (unsigned int i = 0; i < rv.size(); i++){
        rv[i] /= diver;
    }
    return rv;
}

void initRand(std::vector<double> &rv, int size){
    rv.clear();
    rv.reserve(size);
    for (int i = 0; i < size; i++){
        rv.push_back(double(rand()));
    }
}

bool isVectorsNear(std::vector<double> &v1,std::vector<double> &v2){
    for (int i = 0; i < v1.size(); i++){
        if (fabs(v1[i] - v2[i]) > ACCURACY)
            return false;
    }
    return true;
}

void saveVector(std::vector<double> &myVector, std::string fileName) {
    FILE *fd;
    fd = fopen(("result_" + fileName).c_str(), "w+t");
    for (const auto &value : myVector){
        fprintf(fd, "%lf ", value);
    }
    fclose(fd);
}

void sendToChild(std::vector<double> &rv,int size){
    double *rightVector_c = (double*)malloc(sizeof(double) * size);
    for(int i = 0; i < size; i++){
        rightVector_c[i] = rv[i];
    }
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get current process number

    for (int i = 1; i < numProc; i++) {
        // вектор
        MPI_Send((void*)rightVector_c,
                 size,
                 MPI_DOUBLE,
                 i,
                 0,
                 MPI_COMM_WORLD);
//        MPI_Send(&size,
//                 1,
//                 MPI_INT,
//                 i,
//                 1,
//                 MPI_COMM_WORLD);
//        printf("rv is %d",rightVector_c);
    }
    free(rightVector_c);
}

void sendMatrixToChild(Matrix &matr){
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get current process number
    double * strMatrix = (double*)malloc(sizeof(double) * matr.size);
    for (int i = 0; i < matr.size ; i++){
        int dest = (i%(numProc - 1))+1;
        for(int j = 0; j < matr.size; j++){
            strMatrix[j] = matr.values[i * matr.size + j];
        }
//        printf("str is %d",strMatrix);
        // iя - строка матрицы
        MPI_Send(strMatrix,
                 matr.size,
                 MPI_DOUBLE,
                 dest,
                 i/(numProc - 1)+1,
                 MPI_COMM_WORLD);

    }
    free(strMatrix);
}

void collectResults(std::vector<double> &result,int size){
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get current process number

    for (int i = 0; i < size; i++){
        int dest = (i%(numProc - 1))+1;
        double value;
        MPI_Recv(&value, 1, MPI_DOUBLE,
                 dest, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        result[i] = value;
    }
}

int main(int argc, char *argv[])
{


    MPI_Init(&argc,&argv);
    int rank,numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);// get current process number
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get ammount of processes
    if (rank == 0){
        std::string filenameactual = fileNames[CHOICE];
        std::vector<double> myVector;
        printf("start 0 process\n");
        Matrix matr(filenameactual);
        std::vector<double> rv,rvold;
        initRand(rv,matr.size);
        double time1,time2;
        time_t time3,time4;
        time1 = MPI_Wtime();
        time3 = clock();
        for (int i = 0; i < MAX_ITERATION; i++){
                rvold = rv;

            sendToChild(rvold,matr.size);

            sendMatrixToChild(matr);

            collectResults(rv,matr.size);

            rv = _div(rv, _dist(rv));

            if (isVectorsNear(rv,rvold))
                break;
        }
        time2 = MPI_Wtime();
        time4 = clock();
        std::cout << "wtime exec is " << time2 - time1 << std::endl;
        std::cout << "time exec is " << time4 - time3 << std::endl;
        saveVector(rv,filenameactual);

    } else {
        while(1){
            int size = sizeMatrix[CHOICE];

            int countStrs = size / (numProc-1) + (rank <= size % (numProc - 1) ? 1 : 0);
            // 2
            double *rVector = (double*)malloc(sizeof(double)*size);
            double *matrStr = (double*)malloc(sizeof(double)*size);

            MPI_Recv(rVector,size,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            // i - message Tag
            for (int i = 1; i <= countStrs; i++) {

                MPI_Recv(matrStr, size, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                double cellValue = 0;
                // 3 Умножение строки матрицы на вектор
                for (int j = 0; j < size; j++) {
                    cellValue += rVector[j] * matrStr[j];
                }
                // 4
                MPI_Send(&cellValue,
                         1,
                         MPI_DOUBLE,
                         0, // (rank - 1) + (msgTag - 1) * (numProc - 1)
                         (rank - 1) + (i - 1) * (numProc - 1), // снова превратить в номер строки
                         MPI_COMM_WORLD);
            }
        }
    }

    MPI_Finalize();


    return EXIT_SUCCESS;










}