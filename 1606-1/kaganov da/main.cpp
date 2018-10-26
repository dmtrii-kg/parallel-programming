#include <iostream>
#include <random>
#include <ctime>
#include "mpi.h"
using namespace std;

char* FillStr(int n)
{
	int j;
	char* str;
	str = new char[n];
	for(int i = 0; i < n; i++)
		str[i] = rand()%25 + 97;
	for(int i = 0; i < n / 10; i++)
	{
		j = rand() % n;
		str[j] = '.';
	}
	return str;
}

bool searchDot(char c)
{
	switch(c)
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
	int flag;
	int n = 17;						// string lenght
	int counter = 0;					// number of sentences
	int my_counter = 0;
	int remainderDiv = 0;
	int ProcNum, ProcRank, RecvRank;
	double ls_time, le_time, start_time, end_time;

	// random fill of string
	
	// cout << "Введите длину строки: "; cin >> n;
	char* str;
	str = new char[n];
	str = FillStr(n);
	cout << "Полученная строка: " << str << endl;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	// linear block
	if(ProcRank == 0)
	{
		ls_time = MPI_Wtime();
		for(int i = 1; i < n; i++)
		{
			if(searchDot(str[i]))
			{
				if(!searchDot(str[i - 1]))
				{
					counter++;
				}
			}
		}

	// linear result
		le_time = MPI_Wtime();
		cout << "\n LINEAR \n";
		cout << "В введённой строке найдено " << counter << "  предложений\n";
		cout << "За время: " << le_time - ls_time << endl;
	// end linear block
	}

	// parallel block
	start_time = MPI_Wtime();
	
	MPI_Bcast(str, n, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	int size = n / ProcNum;
	int start = size * ProcRank;
	int end = size * (ProcRank + 1);
	
	for(int i = start; i < end; i++)
	{
		if(searchDot(str[i]))
		{
			if(!searchDot(str[i - 1]))
			{
				my_counter++;
			}
		}
	}
	
	MPI_Reduce(&my_counter, &counter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if(ProcRank == 0)
	{
		remainderDiv = n % size;
		if(remainderDiv)
		{
			for(int i = 1; i < remainderDiv; i++)
			{
				if(searchDot(str[n - remainderDiv + i]))
				{
					if(!searchDot(str[n - remainderDiv + i - 1]))
					{
						counter++;
					}
				}
			}
		}
	
	// parallel results
		end_time = MPI_Wtime();
		cout << "\n PARALLEL \n";
		cout << "В введённой строке найдено " << counter << " предложений\n";
		cout << "За время: " << end_time - start_time << endl;;
	}
	MPI_Finalize();
	// end parallel block

	return 0;
}
