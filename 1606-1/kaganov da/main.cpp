#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>
#include "mpi.h"
using namespace std;

void FillStr(char *str, int n)
{
	int j;
	for (int i = 0; i < n; i++)
		str[i] = rand() % 25 + 97;
	str[n] = '\0';
	for (int i = 0; i < n / 10; i++)
	{
		j = 1 + rand() % (n - 1);
		str[j] = '.';
	}
}

bool searchDot(char c)
{
	switch (c)
	{
	case '.':
	case '!':
	case '?':
		return true;
	default:
		return false;
	}
}

int main(int argc, char *argv[])
{
	srand(time(0));
	int n = 1000000;				    // string lenght
	int counter = 0;					// number of sentences
	int my_counter = 0;
	int remainderDiv = 0;
	int ProcNum, ProcRank, RecvRank;
	double ls_time, le_time, start_time, end_time;

	// random fill of string
    int* CounterVec;
    char* str;
    char* my_str;
	str = new char[n + 1];
	FillStr(str, n);

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0)
	{
    // linear block
        // cout << "\nEnter the string: "; cin >> setw(n) >> str;
		cout << "\nReceived string: " << str << endl;

		ls_time = MPI_Wtime();
		for (int i = 1; i < n; i++)
		{
			if (searchDot(str[i]))
			{
				if (!searchDot(str[i - 1]))
				{
					counter++;
				}
			}
		}

    // linear result
		le_time = MPI_Wtime();
		cout << "\n LINEAR \n";
        cout << counter << " sentences found in the entered string\n";
        counter = 0;
		cout << "Time: " << le_time - ls_time << endl;
    // end linear block

    // parallel block
        start_time = MPI_Wtime();
        CounterVec = new int[ProcNum];
        for(int i = 1; i < ProcNum; i++)
            MPI_Send(str + n / ProcNum * i, n / ProcNum, MPI_CHAR, i, 0, MPI_COMM_WORLD);
		
        for (int i = 1; i < n / ProcNum; i++)
        {
            if (searchDot(str[i]))
            {
                if (!searchDot(str[i - 1]))
                {
                    my_counter++;
                }
            }
        }
        CounterVec[0] = my_counter;

        for(int i = 1; i < ProcNum; i++)
            MPI_Recv(&CounterVec[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
    
    // check the remainder of string
        remainderDiv = n % ProcNum;
        if (remainderDiv)
        {
            for (int i = 0; i < remainderDiv; i++)
            {
                if (searchDot(str[n - remainderDiv + i]))
                {
                    if (!searchDot(str[n - remainderDiv + i - 1]))
                    {
                        counter++;
                    }
                }
            }
        }
        CounterVec[0] += counter;

    // parallel results
        counter = 0;
        for(int i = 0; i < ProcNum; i++)
        {
            counter += CounterVec[i];
            // cout << CounterVec[i] << endl;
        }
        delete[]CounterVec;
        end_time = MPI_Wtime();
        cout << "\n PARALLEL \n";
        cout << counter << " sentences found in the entered string\n";
        cout << "Time: " << end_time - start_time << endl << endl;
    } else {
        my_str = new char[n / ProcNum];
        MPI_Recv(my_str, n / ProcNum, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        for (int i = 0; i < n / ProcNum; i++)
        {
            if (searchDot(my_str[i]))
            {
                if (!searchDot(my_str[i - 1]))
                {
                    my_counter++;
                }
            }
        }
        delete[]my_str;
        MPI_Send(&my_counter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
	MPI_Finalize();
	// cin.ignore();
	// end parallel block
    delete[]str;
    
	return 0;
}
