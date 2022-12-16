//
// Created by Dan on 30.11.2022.
//
#include <fcntl.h>
#include "Matrix.h"

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

std::vector<double> Matrix::multiply(std::vector<double> rightVector) {
    int rank,numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);// get current process number
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get ammount of processes
//    printf("size is %d\n",numProc);
    printf("rank is %d of %d\n",rank, numProc);
    std::vector<double> _values(values);
    std::vector<double> res;
    res.resize(size);

    if (0 == rank){
        // 1. Отослать на другие процессы
        // 2. Принять там сообщение
        // 3. Провести работу,
        // 4. Отослать обратно
        // 5. Объединить присланные обратно результаты
        // 6. Вывести

        // 1
        double *rightVector_c = (double*)malloc(sizeof(double) * size);
        for(int i = 0; i < size; i++){
            rightVector_c[i] = rightVector[i];
        }
        for (int i = 1; i < numProc; i++) {
            // вектор
            MPI_Send((void*)rightVector_c,
                    size,
                    MPI_DOUBLE,
                    i,
                    0,
                    MPI_COMM_WORLD);
            printf("rv is %d",rightVector_c);
        }
        free(rightVector_c);

        for (int i = 0; i < size ; i++){
            int dest = (i%(numProc - 1))+1;
//            std::vector<double> strMatrix;
            double * strMatrix = (double*)malloc(sizeof(double) * size);
            for(int j = 0; j < size; j++){
                strMatrix[j] = values[i * size + j];
            }
//            strMatrix.insert(strMatrix.begin(),
//                             values.begin() + i * size,
//                             values.begin() + (i + 1) * size);
            printf("str is %d",strMatrix);
            // iя - строка матрицы
            MPI_Send(strMatrix,
                    size,
                    MPI_DOUBLE,
                    dest,
                    i/(numProc - 1)+1,
                    MPI_COMM_WORLD);
            free(strMatrix);

        }
        // 5
        for (int i = 0; i < size; i++){
            int dest = (i%(numProc - 1))+1;
            double value;
            MPI_Recv(&value, 1, MPI_DOUBLE,
                     dest, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            res[i] = value;
        }
    } else { // Child processes
        int countStrs = size / (numProc-1) + (rank <= size % (numProc - 1) ? 1 : 0);
        // 2
        double *rVector = (double*)malloc(sizeof(double)*size);
        double *matrStr = (double*)malloc(sizeof(double)*size);

        MPI_Recv(rVector,size,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        // i - message Tag
        for (int i = 1; i <= countStrs; i++){

            MPI_Recv(matrStr,size,MPI_DOUBLE,0,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

            double cellValue = 0;
            // 3 Умножение строки матрицы на вектор
            for (int j = 0; j < size; j++){
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
//    MPI_Finalize();
    //6
    return res;
}
