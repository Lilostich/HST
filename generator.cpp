//
// Created by Dan on 30.11.2022.
//
#include <math.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

/// @param[in] fileSize - размер файла в мб
void generateMatrixInFile(std::string fileName, int fileSize){
    int bytesFileSize = fileSize * 1024 * 1024;
    int matrixSize = sqrt(bytesFileSize / 8);
    FILE *fd;

    fd = fopen(fileName.c_str(),"w+b");

    fprintf(fd,"%d ",matrixSize);
    for (int i = 0; i < matrixSize; i++){
        for (int j = 0; j < matrixSize; j++)
            fprintf(fd,"%lf ",double(rand()%1000));
    }
    fclose(fd);
}

int main(int argc,char ** argv) {
    if (argc < 2){
        printf("Usage: generate <size> [size] ...\n");
    } else {
        std::string defaultFileName{"matrixFile"};
        std::string defaultFileEnding{".txt"};
        for (int i = 1; i < argc; i++){
            std::string fileName = defaultFileName + std::string(argv[i]) + defaultFileEnding;
            int size = atoi(argv[i]);
            generateMatrixInFile(fileName, size);
            printf("Success generate file '%s'\n",fileName.c_str());
        }
    }
}