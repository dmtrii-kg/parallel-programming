#include <iostream>
#include <random>
#include <time.h>
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
	int n;							// длина введённой строки
	int counter = 0;					// колличество предложений в строке
	int my_counter = 0;
	int remainderDiv = 0;
	int ProcNum, ProcRank, RecvRank;
	double ls_time, le_time, start_time, end_time;

	// автоматическое заполнение строки
	cout << "Введите длину строки: "; cin >> n;
	char* str;
	str = new char[n];
	str = FillStr(n);
	cout << "Полученная строка: " << str << endl;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	// линейная программа
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
		le_time = MPI_Wtime();
		cout << "В введённой строке найдено " << counter << "  предложений\n";
		cout << "За время: " << le_time - ls_time << endl;
	}

	// параллельная программа
	start_time = MPI_Wtime();
	MPI_Bcast(str, n, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	int size = n / ProcNum;
	int start = size * ProcRank;
	int end = size * (ProcRank + 1);
	if(ProcRank == ProcNum - 1) end = size;
	
	for(int i = start + 1; i < end; i++)
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
		end_time = MPI_Wtime();
		cout << "В введённой строке найдено " << counter << " предложений\n";
		cout << "За время: " << start_time - end_time << endl;;
	}
	
	MPI_Finalize();
	return 0;
}
