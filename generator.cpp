//
// Created by Dan on 30.11.2022.
//
#include <math.h>
#include <fcntl.h>
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
            fprintf(fd,"%lf ",double(rand()));
    }
    fclose(fd);
}

int main(int argc,char ** argv) {
    if (argc != 2){
        printf("Usage: generate <size>\n");
    } else {
        std::string defaultFileName{"matrixFile"};
        std::string defaultFileEnding{".txt"};
        std::string fileName = defaultFileName + std::string(argv[1]) + defaultFileEnding;
        int size = atoi(argv[1]);
        generateMatrixInFile(fileName, size);
        printf("Success generate file '%s'\n",fileName.c_str());
    }
}