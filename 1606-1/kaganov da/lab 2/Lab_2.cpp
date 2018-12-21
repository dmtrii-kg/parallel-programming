#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

void createMatrix(int *matrix, int matrixSize) {
    for (int i = 0; i < matrixSize * matrixSize; i++) {
        matrix[i] = rand() % 10;
    }
}

void printMatrix(int *matrix, int matrixSize) {
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            cout << matrix[i * matrixSize + j] << " ";
        }
        cout << endl;
    }
}

/* void transposition(int *matrix, int *newMatrix, int matrixSize) {
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            newMatrix[i * matrixSize + j] = matrix[j * matrixSize + i];
        }
    }
} */

void linearMultiplyMatrix(int *matrixA, int *matrixB, int *matrixResult, int matrixSize) {
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            matrixResult[i * matrixSize + j] = 0;
            for (int k = 0; k < matrixSize; k++) {
                matrixResult[i * matrixSize + k] += matrixA[i * matrixSize + j] * matrixB[j * matrixSize + k];
            }
        }
    }
}

void parallelMultiplyMatrix(int *pMatrixA, int *pMatrixB, int *pMatrixC, int blockSize, int matrixSize, int tmp, int procRank, int procNum) {
    for (int i = 0; i < blockSize; i++) {
        for (int k = 0; k < matrixSize; k++) {
            for (int j = 0; j < blockSize; j++) {
                pMatrixC[i * matrixSize + j + blockSize * ((tmp + procRank) % procNum)] += pMatrixA[i * matrixSize + k] * pMatrixB[k * blockSize + j];
                //cout << "Матрица А [" << i << "] x [" << j << "] " << pMatrixA[i * matrixSize + k] << endl;
                //cout << "Матрица B [" << i << "] x [" << j << "] " << pMatrixB[k * blockSize + j] << endl;
                //cout << "Матрица С [" << i << "] x [" << j << "] " << pMatrixC[i * matrixSize + j + blockSize * ((tmp + procRank) % procNum)] << endl;
            }
        }
    }
}

void checkResults(int *linearRes, int *parallelRes, int matrixSize) {
    for (int i = 0; i < matrixSize; i++)
        for (int j = 0; j < matrixSize; j++)
            if (linearRes[i * matrixSize + j] != parallelRes[i * matrixSize + j]) {
                cout << "\nError! Linear and parallel results are not equal\n\n";
                return;
            }
    cout << "\nMatrices are equal\n\n";
}

int main(int argc, char *argv[]) {
    srand((int)time(0));
    int matrixSize;
    int blockSize;
    //int partA = 0;
    //int partB = 0;
    //int remainderA = 0;
    //int remainderB = 0;
    int* matrixA = nullptr;
    int* matrixB = nullptr;
    int* matrixC = nullptr;
    //int* lResult = nullptr;
    int* pMatrixA = nullptr;
    int* pMatrixB = nullptr;
    int* pMatrixC = nullptr;
    //int* pResult = nullptr;
    int tmp = 0;
    int procNum = 0;
    int procRank = 0;
    int nextProc = 0;
    int prevProc = 0;
    //int part_aSize = 0;
    //int part_bSize = 0;
    //int part_cSize = 0; 
    //int recvRank = 0;
    //int index = 0;
    double linerStart = 0.0;
    double linerEnd = 0.0;
    double parallelStart = 0.0;
    double parallelEnd = 0.0;
    
    
    // read data from the cmd
    if (argc > 1) {
        matrixSize = atoi(argv[1]);
    } else {
        matrixSize = 4;
    }
    
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    //new data type - columns
    MPI_Datatype MPI_BLOCK;
    
    blockSize = matrixSize / procNum;
    
    if (procRank == 0) {
        // memory allocation for matrices
        matrixA = new int[matrixSize * matrixSize];
        matrixB = new int[matrixSize * matrixSize];
        matrixC = new int[matrixSize * matrixSize];
        // fill the matrices
        createMatrix(matrixA, matrixSize);
        createMatrix(matrixB, matrixSize);
        
        if (matrixSize < 5) {
            cout << "\nMatrix A:\n";
            printMatrix(matrixA, matrixSize);
            cout << "\nMatrix B:\n";
            printMatrix(matrixB, matrixSize);
        }
        
// LINEAR BLOCK
        linerStart = MPI_Wtime();
        linearMultiplyMatrix(matrixA, matrixB, matrixC, matrixSize);
        linerEnd = MPI_Wtime();
        // linear result
        cout << "\n---===LINE===---\n";
        if (matrixSize < 6) {
            cout << "Matrix C:\n";
            printMatrix(matrixC, matrixSize);
        }
        cout << "Time: " << linerEnd - linerStart << " sec\n";
// END LINEAR BLOCK

// PARALLEL BLOCK
        parallelStart = MPI_Wtime();
    }
    pMatrixA = new int[blockSize * matrixSize];
    pMatrixB = new int[blockSize * matrixSize];
    pMatrixC = new int[blockSize * matrixSize];
    for (int i = 0; i < blockSize * matrixSize; i++) {
        pMatrixC[i] = 0;
    }
    
    MPI_Type_vector(matrixSize, blockSize, matrixSize, MPI_INT, &MPI_BLOCK);
    MPI_Type_commit(&MPI_BLOCK);

    if (procRank == 0) {
        for (int i = 1; i < procNum; i++) {
            MPI_Send(matrixA + i * matrixSize * blockSize, matrixSize * blockSize, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(matrixB + i * blockSize, 1, MPI_BLOCK, i, 0, MPI_COMM_WORLD);
        }
        for (int i = 0; i < blockSize * matrixSize; i++) {
            pMatrixA[i] = matrixA[i];
        }
        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < blockSize; j++) {
                pMatrixB[i * blockSize + j] = matrixB[i * matrixSize + j];
            }
        }
    } else {
        MPI_Recv(pMatrixA, blockSize * matrixSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(pMatrixB, 1, MPI_BLOCK, 0, 0, MPI_COMM_WORLD, &status);
        
        /* cout << "Proc = " << procRank << endl;
        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j < blockSize; j++) {
                cout << pMatrixB[i * matrixSize + j] << " ";
            }
            cout << endl;
        } */
    }
    for (tmp = 0; tmp < procNum; tmp++) {
        parallelMultiplyMatrix(pMatrixA, pMatrixB, pMatrixC, blockSize, matrixSize, tmp, procRank, procNum);
        if (tmp < procNum - 1) {
            nextProc = (procRank + 1) % procNum;
            prevProc = (procRank + procNum - 1) % procNum;
            MPI_Sendrecv_replace(pMatrixB, blockSize * matrixSize, MPI_INT, prevProc, 0, nextProc, 0, MPI_COMM_WORLD, &status);
        }
    }
    if (procRank == 0) {
        for (int i = 1; i < procNum; i++) {
            MPI_Recv(matrixC + i * matrixSize * blockSize, matrixSize * blockSize, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        MPI_Send(pMatrixC, blockSize * matrixSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    
    // parallel results
    if (procRank == 0) {
        parallelEnd = MPI_Wtime();
        cout << "\n---===PARALLEL===---\n";
        if (matrixSize < 6) {
            cout << "Matrix C:\n";
            printMatrix(pMatrixC, matrixSize);
        }
        cout << "Time: " << parallelEnd - parallelStart << " sec\n";
        checkResults(matrixC, pMatrixC, matrixSize);
    }

    MPI_Type_free(&MPI_BLOCK);
    MPI_Finalize();
    if (pMatrixA != NULL ) delete []pMatrixA;
    if (pMatrixB != NULL ) delete []pMatrixB;
    if (pMatrixC != NULL ) delete []pMatrixC;
// END PARALLEL BLOCK
    if (matrixA != NULL ) delete []matrixA;
    if (matrixB != NULL ) delete []matrixB;
    if (matrixC != NULL ) delete []matrixC;
    
    return 0;
}
