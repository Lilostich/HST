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

void findMyValue(Matrix &m,std::vector<double> &rv){
    std::vector<double> rvold;
    time_t sum = 0;
    time_t a = clock();
    for (int i = 0; i < MAX_ITERATION; i++){
        rvold = rv;
        rv = _div(m.multiply(rv), _dist(m.multiply(rv)));
        if (isVectorsNear(rv,rvold))
                break;
    }
    time_t b = clock();
    sum += b - a;
    std::cout << "time is " << sum << std::endl;
}

void saveVector(std::vector<double> &myVector, std::string fileName) {
    FILE *fd;
    fd = fopen(("result_" + fileName).c_str(), "w+t");
    for (const auto &value : myVector){
        fprintf(fd, "%lf ", value);
    }
    fclose(fd);
}

int main(int argc, char *argv[])
{
    for (int i = 0; i < FILE_COUNT; i++){
        Matrix matrix(fileNames[i]);
        std::vector<double> myVector;
        initRand(myVector,matrix.getSize());
        findMyValue(matrix,myVector);
        saveVector(myVector,fileNames[i]);
    }
    return EXIT_SUCCESS;
}

// Последовательная версия
/*
 * time is 79
 * time is 333
 * time is 697
 * time is 1580
 * time is 3078
 * */