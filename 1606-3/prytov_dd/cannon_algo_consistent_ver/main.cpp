#include <iostream>
#include <ctime>
#include <omp.h>

using namespace std;

inline int idx(int i, int j, int size) {
	return i * size + j;
}

int matrix_print(double *&M, const int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			cout << M[idx(i, j, size)] << "\t";
		cout << endl;
	}
	cout << endl;
	return 0;
}

int matrix_generate(double *&M, const int size) {
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			M[idx(i, j, size)] = rand() % 100 - 50;
	return 0;
}

int matrix_mult(double *&A, double *&B, double *&C, int size) {
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) {
			C[idx(i, j, size)] = 0;
			for (int k = 0; k < size; k++)
				C[idx(i, j, size)] += A[idx(i, k, size)] * B[idx(k, j, size)];
		}
	return 0;
}

int mult_process_con(double *&A, double *&B, double *&C, int size, int number_blocks_in_row, int i) {
	int block_size = size / number_blocks_in_row;
	int block_row, block_col, row, col, row_begin, row_end, colomn_begin, colomn_end, id_thread;
	double *block_a = nullptr;
	double *block_b = nullptr;
	double *block_c = nullptr;
	for (id_thread = 0; id_thread < number_blocks_in_row * number_blocks_in_row; id_thread++) {
		row_begin = (id_thread / number_blocks_in_row) * block_size;
		row_end = row_begin + block_size;
		colomn_begin = (id_thread % number_blocks_in_row) * block_size;
		colomn_end = colomn_begin + block_size;

		block_a = new double[block_size * block_size];
		block_b = new double[block_size * block_size];
		block_c = new double[block_size * block_size];

		for (row = row_begin, block_row = 0; row < row_end; row++, block_row++) {
			for (col = colomn_begin, block_col = 0; col < colomn_end; col++, block_col++) {
				block_a[idx(block_row, block_col, block_size)] = A[idx(row, (col + i * block_size + (id_thread / number_blocks_in_row) * block_size) % size, size)];
				block_b[idx(block_row, block_col, block_size)] = B[idx((row + i * block_size + (id_thread % number_blocks_in_row) * block_size) % size, col, size)];
			}
		}

		matrix_mult(block_a, block_b, block_c, block_size);

		for (row = row_begin, block_row = 0; row < row_end; row++, block_row++) {
			for (col = colomn_begin, block_col = 0; col < colomn_end; col++, block_col++)
				C[idx(row, col, size)] += block_c[idx(block_row, block_col, block_size)];
		}

		delete[]block_a;
		delete[]block_b;
		delete[]block_c;
	}
	return 0;
}

int cannon_con(double *&A, double *&B, double *&C, int size, int number_blocks_in_row) {
	for (int i = 0; i < number_blocks_in_row; i++) {
		mult_process_con(A, B, C, size, number_blocks_in_row, i);
	}
	return 0;
}

int main(int argc, char** argv) {
	double *A = nullptr;
	double *B = nullptr;
	double *C1 = nullptr;
	double *C2 = nullptr;
	int size = 0;
	int num_thrd = 1;
	int num_blocks_in_row = 1;
	int flag = 1;

	size = atoi(argv[1]);
	if (argc == 3) {
		num_thrd = atoi(argv[2]);
	}
	if (size % int(sqrt(num_thrd)) != 0) {
		cout << "Error! Size must be divided into sqrt(number of threads)" << endl;
		return 0;
	}
	
	num_blocks_in_row = int(sqrt(num_thrd));

	srand((unsigned int)time(NULL));
	A = new double[size*size];
	B = new double[size*size];
	C1 = new double[size*size];
	C2 = new double[size*size];
	for (int i = 0; i < size * size; i++) {
		C1[i] = 0;
		C2[i] = 0;
	}
	matrix_generate(A, size);
	matrix_generate(B, size);

	/*cout << "Matrix A:" << endl;
	matrix_print(A, size);
	cout << endl;
	cout << "Matrix B:" << endl;
	matrix_print(B, size);*/

	const auto start_time_con = omp_get_wtime();
	cannon_con(A, B, C1, size, num_blocks_in_row);
	const auto end_time_con = omp_get_wtime();
	
	/*cout << "Matrix C1:" << endl;
	matrix_print(C1, size);
	cout << endl;*/

	cout << "Cannon consistent time: " << end_time_con - start_time_con << endl;

	return 0;
}