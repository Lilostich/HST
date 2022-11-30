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



Matrix::Matrix(std::vector<double> matr) :values(matr){size = sqrt(matr.size());}

void Matrix::initSize(int n) {values.reserve(n*n); size = n*n;}

double &Matrix::getValue(int i, int j) {return values[j + i * size];}

std::vector<double> Matrix::multiply(std::vector<double> rightVector) {
    std::vector<double> _values(values);
    int _size = this->size;
    if (size != (int)rightVector.size())
        return std::vector<double>();
    else{
        std::vector<double> res;
        res.resize(size);
        int j = 0;
        for (int i = 0; i < _size; i++){
            res[i] = 0;
            for (j = 0; j < _size; j++){
                res[i] += rightVector[j] * _values[i * _size + j];
            }
        }
        return res;
    }
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
