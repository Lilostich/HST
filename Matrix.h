//
// Created by Dan on 30.11.2022.
//

#ifndef MYVECTOR_MATRIX_H
#define MYVECTOR_MATRIX_H

#include <vector>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <fstream>


class Matrix {
public:
    Matrix(std::vector<double> matr);
    Matrix(std::string fileName);
    Matrix(){}
    void initSize(int n);

    double& getValue(int i, int j);

    std::vector<double> multiply(std::vector<double> rightVector);

    void save(std::string fileName);
    int getSize(){return size;}
private:
    std::vector<double> values;
    int size;
};

#endif //MYVECTOR_MATRIX_H
