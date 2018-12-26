//17.Поразрядная сортировка для целых чисел с простым слиянием.
#include <mpi.h>
#include <ctime>
#include <iostream>

#define MainProc 0

using namespace std;

void swp(unsigned int *a, unsigned int *b) {
    unsigned int tmp = *a;
    *a = *b;
    *b = tmp;
}

unsigned int* merge(unsigned int *arr1, unsigned int num1,
                    unsigned int *arr2, unsigned int num2) {
    unsigned int *res;
    unsigned int i = 0,
    j = 0,
    index = 0;
    
    res = new unsigned int[num1 + num2];
    
    while (i < num1 && j < num2) {
        if (arr1[i] < arr2[j])
            res[index++] = arr1[i++];
        else
            res[index++] = arr2[j++];
    }
    
    while (i < num1)
        res[index++] = arr1[i++];
    
    while (j < num2)
        res[index++] = arr2[j++];
    return res;
}

void radixSortMSD(unsigned int *from, unsigned int *to, unsigned int bit) {
    if (!bit || to < from + 1) return;
    
    unsigned int *ll = from, *rr = to - 1;
    
    for (;;) {
        while (ll < rr && !(*ll & bit)) ll++;
        while (ll < rr && (*rr & bit)) rr--;
        if (ll >= rr) break;
        swp(ll, rr);
    }
    
    if (!(bit & *ll) && ll < to) ll++;
    bit >>= 1;
    
    radixSortMSD(from, ll, bit);
    radixSortMSD(ll, to, bit);
}

void radixSortLSD(unsigned int *a, unsigned int count) {
    unsigned int mIndex[4][256] = { 0 };        // count and index matrix
    unsigned int *b = new unsigned int[count];  // allocate temp array
    unsigned int i, j, m, n;
    unsigned int u;
    for (i = 0; i < count; i++) {
        u = a[i];
        for (j = 0; j < 4; j++) {
            mIndex[j][static_cast<unsigned int>(u & 0xff)]++;
            u >>= 8;
        }
    }
    for (j = 0; j < 4; j++) {           // convert to indices
        m = 0;
        for (i = 0; i < 256; i++) {
            n = mIndex[j][i];
            mIndex[j][i] = m;
            m += n;
        }
    }
    for (j = 0; j < 4; j++) {           // radix sort
        for (i = 0; i < count; i++) {   // sort by current lsb
            u = a[i];
            m = static_cast<unsigned int>(u >> (j << 3)) & 0xff;
            b[mIndex[j][m]++] = u;
        }
        swap(a, b);  // swap ptrs
    }
    delete[] b;
}

void Check(unsigned int* arr, unsigned int procNum) {
    for (int i = 0; i < procNum - 1; i++)
        if (arr[i] > arr[i + 1]) {
            cout << "Error! Array not sorted!\n";
            return;
        }
    cout << "Array is sorted\n";
}

void CheckResult(unsigned int* resultL, unsigned int* resultP, unsigned int arrSize) {
    for (int i = 0; i < arrSize; i++)
        if (resultL[i] != resultP[i]) {
            cout << "\nError! Linear and parallel results are not equal\n\n";
            return;
        }
    cout << "\nArrays are the same\n";
}

int main(int argc, char** argv) {
    int status = 0;
    int procRank = 0;
    int procNum = 0;
    int *scounts = NULL;                   // use in scatterv
    int *displs = NULL;                    // use in scatterv
    unsigned int arrSize = 0;
    unsigned int *arr = NULL;              // data
    unsigned int *buff = NULL;             // buffers for message exchanging
    unsigned int *buff2 = NULL;
    unsigned int *resultL = NULL;          // results
    unsigned int *resultP = NULL;
    double linerStart = 0.0;
    double linerEnd = 0.0;
    double parallelStart = 0.0;
    double parallelEnd = 0.0;
    int typeSort = 0;
    int imin = numeric_limits<int>::min();

    status = MPI_Init(&argc, &argv);
    status = MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    status = MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    // read data from the cmd
    if (procRank == MainProc) {
        if (argc > 1) {
            arrSize = atoi(argv[1]);
            if (argc > 2) typeSort = atoi(argv[2]);
        } else {
            arrSize = 10;
            typeSort = 0;
        }
        
        cout << "\nArray size : " << arrSize << endl;
        arr = new unsigned int[arrSize];
        resultL = new unsigned int[arrSize];
        resultP = new unsigned int[arrSize];
        srand((unsigned)time(NULL));
        
        for (int i = 0; i < arrSize; i++)
            arr[i] = resultL[i] = rand() % 200;
        if (arrSize < 31) {
            for (int i = 0; i < arrSize; i++)
                cout << arr[i] << " ";
            cout << endl;
        }
        
// LINEAR BLOCK
        linerStart = MPI_Wtime();
        if (typeSort == 0)
            radixSortMSD(resultL, resultL + arrSize, imin);
        else
            radixSortLSD(resultL, arrSize);
        linerEnd = MPI_Wtime();
        if (arrSize < 31) {
            cout << endl;
            for (int i = 0; i < arrSize; i++)
                cout << resultL[i] << " ";
        }
        cout << "\nLiner time: " << linerEnd - linerStart << endl;
        Check(resultL, arrSize);
// END LINEAR BLOCK
        
// PARALLEL BLOCK
        buff = new unsigned int[arrSize / procNum + arrSize % procNum];
        buff2 = new unsigned int[arrSize / procNum];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    parallelStart = MPI_Wtime();
    MPI_Bcast(&arrSize, 1, MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);
    MPI_Bcast(&typeSort, 1, MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);
    
    scounts = new int[procNum];
    displs = new int[procNum];
    displs[0] = 0;
    scounts[0] = arrSize / procNum + arrSize % procNum;
    
    for (int i = 1; i < procNum; i++) {
        scounts[i] = arrSize / procNum;
        displs[i] = scounts[0] + (i-1)*scounts[i];
    }
    
    if (procRank != 0)
        buff = new unsigned int[scounts[procRank]];
    
    MPI_Scatterv(arr, scounts, displs, MPI_UNSIGNED, buff, scounts[procRank], MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);
    
    if (typeSort == 0)
        radixSortMSD(buff, buff + scounts[procRank], imin);
    else
        radixSortLSD(buff, scounts[procRank]);
    
    if (procRank != 0)
        MPI_Send(buff, scounts[procRank], MPI_UNSIGNED, MainProc, procRank, MPI_COMM_WORLD);
    
    
    if (procRank == 0) {
        for (int i = 0; i < static_cast<unsigned int>(scounts[0]); i++)
            resultP[i] = buff[i];
        for (int i = 1; i < static_cast<unsigned int>(procNum); i++) {
            MPI_Recv(buff2, arrSize / procNum, MPI_UNSIGNED, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            resultP = merge(buff2, arrSize / procNum, resultP, scounts[0] + (i - 1)*scounts[1]);
        }
        // parallel results
        parallelEnd = MPI_Wtime();
        if (arrSize < 31) {
            cout << endl;
            for (unsigned int i = 0; i < arrSize; i++)
                cout << resultP[i] << " ";
        }
        cout << "\nParallel Time: " << parallelEnd - parallelStart << endl;
        Check(resultP, arrSize);
        CheckResult(resultL, resultP, arrSize);
        cout << "\nAverage acceleration: " << (linerEnd - linerStart) / (parallelEnd - parallelStart) << endl;
        cout << endl;
        
        delete[] buff2;
        delete[] resultL;
        delete[] resultP;
        delete[] arr;
    }
    delete[] buff;
    delete[] scounts;
    delete[] displs;
    status = MPI_Finalize();
    
    return 0;
}
