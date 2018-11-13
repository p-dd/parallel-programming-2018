#include <iostream>
#include "mpi.h"
#include <chrono>
#include <random>
#include <string>
using namespace std;

void printMatrix(int size,double* matr)
{
	cout << "Generated matrix:" << endl << endl;
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
			cout << matr[size*i + j] << " | ";
		cout << endl << endl;
	}
}

int main(int argc, char* argv[])
{
	int rank, proc;
	double time = 0.0;
	double* matrix = nullptr;
	double* ans = nullptr;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc);

	const int size = stoi(string(argv[1]));

	bool showMatr = false;
	if (argc == 3 && string(argv[2]) == "show")
		showMatr = true;

	if (rank == 0)
	{
		default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());

		const uniform_real_distribution <double> distribution(-100, 100);

		matrix = new double[size*size];

		cout.precision(5);
		for (int i = 0; i < size*size; ++i)
			matrix[i] = distribution(generator);

		if (showMatr)
		{
			cout << endl;
			printMatrix(size, matrix);
		}

		time = MPI_Wtime();
	}
	if (proc > 1)
	{
		const int nRows = size / proc;
		int residue = size % proc;
		double* mas = residue == 0 ? new double[size*nRows] : new double[size*nRows + size];

		int end = residue == 0 ? size * nRows : size * nRows + size;

		for (int i = 0; i < end; ++i)
		{
			mas[i] = 0;
		}

		int sumEl = 0;
		auto sendcounts = new int[proc];
		auto displs = new int[proc];

		for (int i = 0; i < proc; i++)
		{
			sendcounts[i] = nRows * size;
			if (residue > 0) {
				sendcounts[i] += size;
				residue--;
			}
			displs[i] = sumEl;
			sumEl += sendcounts[i];
		}
		const int myNum = sendcounts[rank];

		MPI_Scatterv(matrix, sendcounts, displs, MPI_DOUBLE, mas, myNum, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		auto res = size % proc == 0 ? new double[nRows] : new double[nRows + 1];

		const int count = size % proc == 0 ? nRows : nRows + 1;

		for (int i = 0; i < count; ++i)
		{
			double s = 0.0;
			for (int j = 0; j < size; ++j)
			{
				s += mas[i*size + j];
			}
			res[i] = s;
		}

		if (rank == 0)
		{
			ans = new double[size];
		}

		for (int i = 0; i < proc; i++)
		{
			sendcounts[i] /= size;
			displs[i] /= size;
		}

		int num = 0;
		while (res[num] != 0 && num != count)
		{
			++num;
		}

		MPI_Gatherv(res, num, MPI_DOUBLE, ans, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0)
		{
			time = MPI_Wtime() - time;
			if (showMatr)
			{
				cout << endl << "Sum of elems in each row:" << endl;
				for (int i = 0; i < size; ++i)
				{
					cout << ans[i] << " | ";
				}
				cout << endl;
			}

			cout << endl << "Parallel algo time: " << time << "sec" << endl;
		}
	}
	else
	{
		if (rank == 0)
		{
			time = MPI_Wtime();
			double* ansIter = new double[size];
			for (int i = 0; i < size; ++i) {
				double s = 0;
				for (int j = 0; j < size; ++j) {
					s += matrix[i*size + j];
				}
				ansIter[i] = s;
			}
			if (showMatr) {
				cout << endl << "Sum of elems in each row:" << endl;
				for (int i = 0; i < size; ++i) {
					cout << ansIter[i] << " | ";
				}
				cout << endl;
			}
			cout << endl << "Iterative algo time: " << MPI_Wtime() - time << "sec" << endl << endl;
		}
	}
	MPI_Finalize();
	return 0;
}
