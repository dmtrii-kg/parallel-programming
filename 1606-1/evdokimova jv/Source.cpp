#include <iostream>
#include <mpi.h>
#include <string>
using namespace std;

void printMatrix(int *data, int row, int col) {
	if (row < 15 && col < 15) {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				cout << data[i * col + j] << " ";
			}
			cout << endl;
		}
	}
}

int* createMatrix(int row, int col) {
	int *matrix;
	matrix = new int[row*col];
	return matrix;
}

void fullMatrix(int *matrix, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			matrix[i*col + j] = rand() % 100;
		}
	}
}
int minSearch(int a, int b) {
	if (a<=b)
		return a;
	else return b;
}
int main(int argc, char **argv) {
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	const int rows = stoi(string(argv[1]));
	const int cols = stoi(string(argv[2]));

	int *matrix = nullptr;
	if (rank == 0) {

		matrix = new int[rows * cols];
		fullMatrix(matrix, rows, cols);
		printMatrix(matrix, rows, cols);
	}

	int partSize = rows / size;
	int *vec = new int[cols * partSize];
	MPI_Scatter(matrix, cols * partSize, MPI_INT, vec, cols*partSize, MPI_INT, 0, MPI_COMM_WORLD);

	int *localMin = new int[partSize];

	for (int i = 0; i < partSize; i++) {
		int min = INT_MAX;
		for (int j = 0; j < cols; j++) {
			if (vec[i * cols + j] < min) {
				min= vec[i * cols + j];
			}
		}
		localMin[i] = min;
	}
	int *totalMin = nullptr;
	if (rank == 0) {
		totalMin = new int[rows];
	}
	MPI_Gather(localMin, partSize, MPI_INT, totalMin, partSize, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		int tail = rows - size * partSize;
		for (int i = tail + 1; i < rows; i++) {
			int min = INT_MAX;
			for (int j = 0; j < cols; j++) {
				min = minSearch(matrix[i*cols + j], min);
			}
			totalMin[i] = min;
		}

		for (int i = 0; i < rows; i++) {
			cout << totalMin[i] << endl;
		}
	}

	MPI_Finalize();
	return 0;
}