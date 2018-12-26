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

/*
void transposition(int *matrix, int *newMatrix, int rows, int columns) {
    for (int i = 0; i < columns; i++) {
        for (int j = 0; j < rows; j++) {
            newMatrix[i * rows + j] = matrix[j * columns + i];
        }
    }
}
*/

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

int* generateOffsetArray(const int columns, const int rows, int numProc, int* amountArray) {
    int* arr = new int[numProc];
    int currOffset = 0;
    arr[0] = 0;
    for (int i = 1; i < numProc; i++) {
        currOffset = arr[i - 1];
        arr[i] = currOffset + amountArray[i - 1];
    }
    return arr;
}

int* generateAmountArray(const int columns, const int rows, int numProc) {
    const int size = columns * rows;
    int curr = 0;
    int i = 0;
    int* arr = new int[numProc];
    for (int j = 0; j < numProc; j++)
        arr[j] = 0;
    while (curr < size) {
        arr[i++ % numProc]++;
        curr += rows;
    }
    return arr;
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
    int* pMatrixRes = nullptr;
    int* pMatrixC = nullptr;
    int tmp = 0;
    int aSize = aRows * aColumns;
    int bSize = bRows * bColumns;
    int part_aSize = 0;
    int part_bSize = 0;
    int part_cSize = 0;
    int procNum = 0;
    int procRank = 0;
    int nextProc = 0;
    int prevProc = 0;
    int index = 0;                      // смещение по столбцу в ленте относительно пришедшего столбца
    int *amountAray = nullptr;
    int *offsetAray = nullptr;
    int *elementsPerProcGather = nullptr;
    int *shiftsGather = nullptr;
    int *recvbuf = nullptr;
    int *resultBuffer = nullptr;
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
    
    MPI_Datatype column;
    MPI_Datatype columnType;
    
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
        // transposition(matrixB, tMatrixB, bRows, bColumns);
        for (int i = aRows * aColumns; i < partA * procNum * aColumns; i++)
            matrixA[i] = 0;
        for (int i = bRows * bColumns; i < partB * procNum * bRows; i++)
            tMatrixB[i] = 0;
        
        if ((aColumns < 5) && (aRows < 5) && (bColumns < 5) && (bRows < 5)) {
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
    pMatrixRes = new int[part_bSize];
    for (int i = 0; i < part_aSize; i++)
        pMatrixA[i] = 0;
    for (int i = 0; i < part_bSize; i++)
        pMatrixB[i] = 0;
    for (int i = 0; i < part_cSize; i++)
        pMatrixC[i] = 0;
    for (int i = 0; i < part_bSize; i++)
        pMatrixRes[i] = 0;
    
    amountAray = new int[procNum];
    offsetAray = new int[procNum];
    
    if (procRank == 0) {
        amountAray = generateAmountArray(bColumns, bRows, procNum);
        offsetAray = generateOffsetArray(bColumns, bRows, procNum, amountAray);
        recvbuf = new int[amountAray[0] * bRows];
        for (int j = 0; j < amountAray[0] * bRows; j++)
            recvbuf[j] = 0;
        //amountAray[0] += remainderB * bColumns;
    }
    
    MPI_Type_vector(bRows, 1, bColumns, MPI_INT, &column);
    MPI_Type_commit(&column);
    
    MPI_Type_create_resized(column, 0, 1 * sizeof(int), &columnType);
    MPI_Type_commit(&columnType);
    
    MPI_Scatter(matrixA, part_aSize, MPI_INT, pMatrixA, part_aSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(matrixB, amountAray, offsetAray, columnType, pMatrixB, part_bSize, columnType, 0, MPI_COMM_WORLD);
    
    /*
    //cout parallel matrixB
    cout << "\nprocNum: " << procRank << endl;
    for (int i = 0; i < bRows; i++) {
        for (int j = 0; j < partB; j++) {
            cout << pMatrixB[i * bRows + j] << " ";
        }
        cout << endl;
    }

    cout << endl;
    for (int i = 0; i < part_bSize * 2; i++) {
        cout << pMatrixB[i] << " ";
    }
    cout << endl;
    */
    
    //create matrix Result
    for (int i = 0, k = 0; i <= 2 * part_bSize / partB - 2; i += 2) {
        for (int j = 0; j < partB; j++) {
            pMatrixRes[k++] = pMatrixB[i * partB + j];
        }
    }
    
    /*
    cout << endl;
    for (int i = 0; i < part_bSize * 2; i++) {
        cout << pMatrixRes[i] << " ";
    }
    cout << endl;
    */
    
    // calculation of elements on the main diagonal
    for (int i = 0; i < partA; i++) {
        for (int j = 0; j < partB; j++) {
            pMatrixC[i * partB * procNum + j + partB * procRank] = 0;
            for (int k = 0; k < aColumns; k++) {
                pMatrixC[i * partB * procNum + j + partB * procRank] += pMatrixA[i * aColumns + k] * pMatrixRes[j + partB * k];
            }
            //cout << "pMatrixC = " << pMatrixC[i * partB * procNum + j + partB * procRank] << endl;
        }
    }
    
    
    nextProc = procRank + 1;
    if (procRank == procNum - 1)
        nextProc = 0;
    prevProc = procRank - 1;
    if (procRank == 0)
        prevProc = procNum - 1;
    
    // cyclic exchange of columns matrix B
    
    for (int i = 0; i < part_bSize; i++) {
        pMatrixRes[i] = 0;
    }
    for (int num = 1; num < procNum; num++) {                                    // num -- the number of perfect shipments
        
        MPI_Sendrecv_replace(pMatrixB, part_bSize, columnType, nextProc, 0, prevProc, 0, MPI_COMM_WORLD, &status);
        
        cout << "\nprocNum: " << procRank << endl;
        for (int i = 0; i < bRows; i++) {
            for (int j = 0; j < partB; j++) {
                cout << pMatrixB[i * bRows + j] << " ";
            }
            cout << endl;
        }
        
        for (int i = 0; i < partA; i++) {
            for (int j = 0; j < partB; j++) {
                tmp = 0;
                for (int k = 0; k < aColumns; k++) {
                    tmp += pMatrixA[i * aColumns + k] * pMatrixRes[j + bRows * k];
                }
                if (procRank - num >= 0) {
                    index = procRank - num;
                } else {
                    index = (procRank - num + procNum);
                }
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
    if (pMatrixA != NULL) delete []pMatrixA;
    if (pMatrixB != NULL) delete []pMatrixB;
    if (pMatrixC != NULL) delete []pMatrixC;
    if (pMatrixRes != NULL) delete []pMatrixRes;
    // END PARALLEL BLOCK
    if (matrixA != NULL) delete []matrixA;
    if (matrixB != NULL) delete []matrixB;
    if (lMatrixC != NULL) delete []lMatrixC;
    if (matrixC != NULL) delete []matrixC;
    if (tMatrixB != NULL) delete []tMatrixB;
    
    return 0;
}
