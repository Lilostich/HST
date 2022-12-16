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
    Matrix(std::string fileName);
    Matrix(){}


    std::vector<double> multiply(std::vector<double> rightVector);

    std::vector<double> values;
    int size;
private:
};

#endif //MYVECTOR_MATRIX_H
