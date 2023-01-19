#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <fstream>


class Matrix {
public:
    Matrix(std::string fileName){
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
    Matrix(){}

    std::vector<double> multiply(std::vector<double> rightVector){
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
    //6
    return res;
}

    std::vector<double> values;
    int size;
};


extern double multiply_str_str(double *leftVector, double *rightVector,int size, time_t *time);

const double ACCURACY = 0.0000000000001;
const std::string fileNames[] {"matrixFile10.txt",
                         "matrixFile50.txt",
                         "matrixFile100.txt",
                         "matrixFile250.txt",
                         "matrixFile500.txt"};
const size_t FILE_COUNT = 5;
const size_t MAX_ITERATION = 1000;

const int sizeMatrix[] {1144,2560,3620, 5724, 8095};

#define CHOICE 0

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
    int a = 0;
    // printf("%d\n", a++);
    for(int i = 0; i < size; i++){
        rightVector_c[i] = rv[i];
    }
    // printf("%d\n", a++);
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get current process number
    // printf("%d\t%d\n",a++, numProc);

    for (int i = 1; i < numProc; i++) {
        // вектор
        // printf("\t%d\n",i);
        MPI_Send((void*)rightVector_c,
                 size,
                 MPI_DOUBLE,
                 i,
                 0,
                 MPI_COMM_WORLD);
        // printf("Send right vector %d->%d as [%lf, %lf,...](size is %d)\n",
                // 0,i,rightVector_c[0],rightVector_c[1],size);
    }
    // printf("%d\n", a);

    free(rightVector_c);
}

void sendMatrixToChild(Matrix &matr){
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get current process number
    double * strMatrix = (double*)malloc(sizeof(double) * matr.size);
    for (int i = 0; i < matr.size ; i++){
        int dest = (i % (numProc - 1)) + 1;
        for(int j = 0; j < matr.size; j++){
            strMatrix[j] = matr.values[i * matr.size + j];
        }
        // iя - строка матрицы
        MPI_Send(strMatrix,
                 matr.size,
                 MPI_DOUBLE,
                 dest,
                 i / (numProc - 1) + 1,
                 MPI_COMM_WORLD);
        // if (i < 10)
            // printf("Send left vector %d->%d as [%lf, %lf,...](size is %d)\n",
                // 0, dest, strMatrix[0], strMatrix[1], matr.size);

    }
    free(strMatrix);
}

void collectResults(std::vector<double> &result,int size){
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get current process number

    for (int i = 0; i < size; i++){
        int dest = (i%(numProc - 1)) + 1;
        double value;
        MPI_Recv(&value, 1, MPI_DOUBLE,
                 dest, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (i < 10)
        // printf("Collect result from %d(%d) proc (value is %lf)\n",dest,i,value);
        result[i] = value;
    }
    // printf("all collected!\n");
}

int main(int argc, char *argv[])
{

    MPI_Init(&argc,&argv);
    int rank,numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // get current process number
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);// get ammount of processes
    printf("rank is %d, numProc is %d\n", rank, numProc);

    if (rank == 0){
        
        std::string filenameactual;
        if (argc > 1)
            filenameactual = argv[1];
        else
            filenameactual = fileNames[CHOICE];
        std::vector<double> myVector;
        // printf("start 0 process\n");
        Matrix matr(filenameactual);
        std::vector<double> rv,rvold;
        initRand(rv,matr.size);
        double time1,time2;
        time_t time3,time4;
        time1 = MPI_Wtime();
        time3 = clock();
        bool first_time = true;
        for (int i = 0; i < MAX_ITERATION; i++){
            // printf("iteration %d\n",i);    
            rvold = rv;
            // printf("rvold = rv\n",i);    

            if (first_time){
                first_time = false;
                
            }
            sendToChild(rvold,matr.size);
            // printf("sendTochild\n",i);    

            sendMatrixToChild(matr);
            // printf("sendMatrixTochild\n",i);    

            collectResults(rv,matr.size);
            // printf("Collect success\n",i);    

            double dist = _dist(rv);
            // printf("dist is %lf\n",dist);

            rv = _div(rv, dist);

            // printf("Vector rv after div %d iteration is [%lf, %lf, ...]]", i , rv[0], rv[1]);

            if (isVectorsNear(rv,rvold)){
                // printf("near!\n");
                break;
            }
            if (i == 999)
                printf("999!\n");
        }
        time2 = MPI_Wtime();
        time4 = clock();
        std::cout << "wtime exec is " << time2 - time1 << std::endl;
        std::cout << "time exec is " << time4 - time3 << std::endl;
        saveVector(rv,filenameactual);

    } else {
        
        while(1){
            bool onlyFirst = true;
            int st = 0;
            int size = sizeMatrix[CHOICE];

            if (argc > 2)
                size = atoi(argv[2]);
            else
                size = sizeMatrix[CHOICE];

            int countStrs = size / (numProc-1) + (rank <= size % (numProc - 1) ? 1 : 0);
            // 2
            // printf("Proc %d step %d\n",rank,st++);

            double *rVector = (double*)malloc(sizeof(double)*size);
            double *matrStr = (double*)malloc(sizeof(double)*size);
            // printf("Proc %d step %d\n",rank,st++);

            MPI_Recv(rVector,size,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            // printf("Proc %d take rVector as [%lf, %lf, ...]\n", rank, rVector[0], rVector[1]);

            // i - message Tag
            // printf("Proc %d step %d\n",rank,st++);
            for (int i = 1; i <= countStrs; i++) {

                MPI_Recv(matrStr, size, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // if (i < 11)
                    // printf("Proc %d take lvector as [%lf, %lf, ...]\n", rank, matrStr[0], matrStr[1]);

                double cellValue = 0;
                // 3 Умножение строки матрицы на вектор
                time_t time;
                cellValue = multiply_str_str(rVector, matrStr,size,&time);
//                for (int j = 0; j < size; j++) {
//                    cellValue += rVector[j] * matrStr[j];
//                }
                if (onlyFirst) {
                    // printf("Proc %d solve cell value as %lf\n", rank, cellValue);
                    onlyFirst = false;
                }
                // 4
                MPI_Send(&cellValue,
                         1,
                         MPI_DOUBLE,
                         0, // (rank - 1) + (msgTag - 1) * (numProc - 1)
                         (rank - 1) + (i - 1) * (numProc - 1), // снова превратить в номер строки
                         MPI_COMM_WORLD);
            
            }
            // printf("Proc %d last step %d\n",rank,st++);
        }
    }

    MPI_Finalize();
    return 0;
}