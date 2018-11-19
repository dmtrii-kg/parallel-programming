#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

void createMatrix(int *matrix, int rows, int columns) {
    for (int i = 0; i < rows * columns; i++) {
        matrix[i] = rand() % 10;
    }
}

void printMatrix(int *matrix, int rows, int columns, int remainderA, int remainderB) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            cout << matrix[i * columns + j + i * remainderB] << " ";
        }
        cout << endl;
    }
}


void transposition(int *matrix, int *newMatrix, int rows, int columns) {
    for (int i = 0; i < columns; i++) {
        for (int j = 0; j < rows; j++) {
            newMatrix[i * rows + j] = matrix[j * columns + i];
        }
    }
}

void multiplyMatrix(int *matrixA, int *matrixB, int *matrixResult, int aRows, int aColumns, int bRows, int bColumns) {
    for (int i = 0; i < aRows; i++) {
        for (int j = 0; j < bColumns; j++) {
            matrixResult[i * bColumns + j] = 0;
            for (int k = 0; k < aColumns; k++) {
                matrixResult[i * bColumns + j] += matrixA[i * aColumns + k] * matrixB[k * bColumns + j];
            }
        }
    }
}

void checkResults(int *linearRes, int *parallelRes, int rows, int columns, int remainderB) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
            if (linearRes[i * columns + j] != parallelRes[i * columns + j + i * remainderB]) {
                cout << "\nError! Linear and parallel results are not equal\n\n";
                return;
            }
    cout << "\nMatrices are equal\n\n";
}

int main(int argc, char *argv[]) {
    srand((int)time(0));
    int aRows = 0;
    int bRows = 0;
    int aColumns = 0;
    int bColumns = 0;
    int partA = 0;
    int partB = 0;
    int remainderA = 0;
    int remainderB = 0;
    int* matrixA = nullptr;
    int* matrixB = nullptr;
    int* matrixC = nullptr;
    int* tMatrixB = nullptr;
    int* lMatrixC = nullptr;
    int* pMatrixA = nullptr;
    int* pMatrixB = nullptr;
    int* pMatrixC = nullptr;
    int tmp = 0;
    // int aSize = aRows * aColumns;
    // int bSize = bRows * bColumns;
    int part_aSize = 0;
    int part_bSize = 0;
    int part_cSize = 0;
    int procNum = 0;
    int procRank = 0;
    int recvRank = 0;
    int nextProc = 0;
    int prevProc = 0;
    int index = 0;                      // смещение по столбцу в ленте относительно пришедшего столбца
    double linerStart = 0.0;
    double linerEnd = 0.0;
    double parallelStart = 0.0;
    double parallelEnd = 0.0;
    
    
    // read data from the cmd
    if (argc > 4) {
        aRows = atoi(argv[1]);
        aColumns = atoi(argv[2]);
        bRows = atoi(argv[3]);
        bColumns = atoi(argv[4]);
    } else {
        aRows = 2;
        aColumns = 3;
        bRows = 3;
        bColumns = 2;
    }
    
    // check matrix sizes
    if (aColumns != bRows) {
        cout << "ERROR::Incorrect matrix values\n";
        return -1;
    }
    
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    partA = static_cast<int>(std::ceil(static_cast<double>(aRows) / static_cast<double>(procNum)));
    partB = static_cast<int>(std::ceil(static_cast<double>(bColumns) / static_cast<double>(procNum)));
    part_aSize = partA * aColumns;
    part_bSize = partB * bRows;
    part_cSize = partA * partB * procNum;
    remainderA = partA * procNum - aRows;
    remainderB = partB * procNum - bColumns;
    
    if (procRank == 0) {
        // memory allocation for matrices
        matrixA = new int[partA * procNum * aColumns];
        matrixB = new int[partB * procNum * bRows];
        matrixC = new int[partA * procNum * partB * procNum];
        tMatrixB = new int[bRows * partB * procNum];
        lMatrixC = new int[aRows * bColumns];
        // fill the matrices and transposition matrix B
        createMatrix(matrixA, aRows, aColumns);
        createMatrix(matrixB, bRows, bColumns);
        transposition(matrixB, tMatrixB, bRows, bColumns);
        for (int i = aRows * aColumns; i < partA * procNum * aColumns; i++)
            matrixA[i] = 0;
        for (int i = bRows * bColumns; i < partB * procNum * bRows; i++)
            tMatrixB[i] = 0;
        
        if ((aColumns < 4) && (aRows < 4) && (bColumns < 4) && (bRows < 4)) {
            cout << "\nMatrix A:\n";
            printMatrix(matrixA, aRows, aColumns, 0, 0);
            cout << "\nMatrix B:\n";
            printMatrix(matrixB, bRows, bColumns, 0, 0);
        }
        
// LINEAR BLOCK
        linerStart = MPI_Wtime();
        multiplyMatrix(matrixA, matrixB, lMatrixC, aRows, aColumns, bRows, bColumns);
        linerEnd = MPI_Wtime();
        // linear result
        cout << "\n---===LINE===---\n";
        if ((aRows < 5) && (bColumns < 5)) {
            cout << "Matrix C:\n";
            printMatrix(lMatrixC, aRows, bColumns, 0, 0);
        }
        cout << "Time: " << linerEnd - linerStart << " sec\n";
// END LINEAR BLOCK

// PARALLEL BLOCK
        parallelStart = MPI_Wtime();
    }
    pMatrixA = new int[part_aSize];
    pMatrixB = new int[part_bSize];
    pMatrixC = new int[part_cSize];
    for (int i = 0; i < part_aSize; i++)
        pMatrixA[i] = 0;
    for (int i = 0; i < part_bSize; i++)
        pMatrixB[i] = 0;
    for (int i = 0; i < part_cSize; i++)
        pMatrixC[i] = 0;
    
    MPI_Scatter(matrixA, part_aSize, MPI_INT, pMatrixA, part_aSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(tMatrixB, part_bSize, MPI_INT, pMatrixB, part_bSize, MPI_INT, 0, MPI_COMM_WORLD);
    // calculation of elements on the main diagonal
    for (int i = 0; i < partA; i++) {
        for (int j = 0; j < partB; j++) {
            pMatrixC[i * partB * procNum + j + partB * procRank] = 0;
            for (int k = 0; k < aColumns; k++) {
                pMatrixC[i * partB * procNum + j + partB * procRank] += pMatrixA[i * aColumns + k] * pMatrixB[j * bRows + k];
            }
        }
    }
    nextProc = procRank + 1;
    if (procRank == procNum - 1)
        nextProc = 0;
    prevProc = procRank - 1;
    if (procRank == 0)
        prevProc = procNum - 1;
    
    // cyclic exchange of columns matrix B
    for (int num = 1; num < procNum; num++) {                                    // num -- the number of perfect shipments
        MPI_Sendrecv_replace(pMatrixB, part_bSize, MPI_INT, nextProc, 0,
                             prevProc, 0, MPI_COMM_WORLD, &status);
        for (int i = 0; i < partA; i++) {
            for (int j = 0; j < partB; j++) {
                tmp = 0;
                for (int k = 0; k < aColumns; k++)
                    tmp += pMatrixA[i * aColumns + k] * pMatrixB[j * bRows + k];
                if (procRank - num >= 0)
                    index = procRank - num;
                else
                    index = (procRank - num + procNum);
                pMatrixC[i * partB * procNum + j + index * partB] = tmp;
            }
        }
    }
    
    // assembly of the resulting matrix
    MPI_Gather(pMatrixC, part_cSize, MPI_INT, matrixC, part_cSize,
               MPI_INT, 0, MPI_COMM_WORLD);
    
    // parallel results
    if (procRank == 0) {
        parallelEnd = MPI_Wtime();
        cout << "\n---===PARALLEL===---\n";
        if ((aRows < 5) && (bColumns < 5)) {
            cout << "Matrix C:\n";
            printMatrix(matrixC, aRows, bColumns, remainderA, remainderB);
        }
        cout << "Time: " << parallelEnd - parallelStart << " sec\n";
        checkResults(lMatrixC, matrixC, aRows, bColumns, remainderB);
    }
    MPI_Finalize();
    if (pMatrixA != NULL ) delete []pMatrixA;
    if (pMatrixB != NULL ) delete []pMatrixB;
    if (pMatrixC != NULL ) delete []pMatrixC;
    // END PARALLEL BLOCK
    if (matrixA != NULL ) delete []matrixA;
    if (matrixB != NULL ) delete []matrixB;
    if (lMatrixC != NULL ) delete []lMatrixC;
    if (matrixC != NULL ) delete []matrixC;
    if (tMatrixB != NULL ) delete []tMatrixB;
    
    return 0;
}
