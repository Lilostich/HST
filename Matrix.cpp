//
// Created by Dan on 30.11.2022.
//
#include <fcntl.h>
#include <mpi.h>

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



Matrix::Matrix(std::vector<double> matr) :values(matr){size = sqrt(matr.size());}

void Matrix::initSize(int n) {values.reserve(n*n); size = n*n;}

double &Matrix::getValue(int i, int j) {return values[j + i * size];}

std::vector<double> Matrix::multiply(std::vector<double> rightVector) {
    int rank,numProc;
    MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);// get current process number
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get ammount of processes
    std::vector<double> _values(values);
    int _size = this->size;
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
        for (int i = 1; i < numProc; i++) {
            // вектор
            MPI_Send(rightVector.data(),
                    size,
                    MPI_DOUBLE,
                    i,
                    0,
                    MPI_COMM_WORLD);
        }
        for (int i = 0; i < size ; i++){
            int dest = (i%(numProc - 1))+1;
            std::vector<double> strMatrix;
            strMatrix.insert(strMatrix.begin(),
                             values.begin() + i * size,
                             values.begin() + (i + 1) * size);
            // iя - строка матрицы
            MPI_Send(strMatrix.data(),
                    size,
                    MPI_DOUBLE,
                    dest,
                    i/(numProc - 1)+1,
                    MPI_COMM_WORLD);

        }
        // 5
        for (int i = 0; i < size; i++){
            int dest = (i%(numProc - 1))+1;
            MPI_Recv(res[i], 1, MPI_DOUBLE,
                     dest, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
            MPI_Send(cellValue,
                    1,
                    MPI_DOUBLE,
                    0, // (rank - 1) + (msgTag - 1) * (numProc - 1)
                    (rank - 1) + (i - 1) * (numProc - 1), // снова превратить в номер строки
                    MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
    //6
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
