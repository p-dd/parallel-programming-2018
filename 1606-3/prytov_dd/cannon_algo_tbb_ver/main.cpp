#include <iostream>
#include <ctime>
#include <omp.h>
#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"

using namespace std;
using namespace tbb;

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

int mult_process_par(double *&A, double *&B, double *&C, int size, int number_blocks_in_row, int i) {
	int block_size = size / number_blocks_in_row;
	int block_row, block_col, row, col, row_begin, row_end, colomn_begin, colomn_end, id_thread;
	double *block_a = nullptr;
	double *block_b = nullptr;
	double *block_c = nullptr;
#pragma omp parallel default(none) private(block_row, block_col, row, col, row_begin, row_end, colomn_begin, colomn_end, id_thread, block_a, block_b, block_c) shared(A, B, C, size, block_size, number_blocks_in_row, i) num_threads(number_blocks_in_row * number_blocks_in_row)
	{
		id_thread = omp_get_thread_num();
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

int cannon_par(double *&A, double *&B, double *&C, int size, int number_blocks_in_row) {
	for (int i = 0; i < number_blocks_in_row; i++) {
		mult_process_par(A, B, C, size, number_blocks_in_row, i);
	}
	return 0;
}

class CannonTask {
private:
	int row_start_idx;
	int col_start_idx;
	int size;
	double *A;
	double *B;
	double *C;
	int block_size;
	int number_blocks_in_row;
public:
	CannonTask(double *_A, double *_B, double *_C,
		int _row_start_idx, int _col_start_idx, int _size, int _block_size, int _number_blocks_in_row) :
		row_start_idx(_row_start_idx), col_start_idx(_col_start_idx), size(_size), A(_A), B(_B),
		C(_C), block_size(_block_size), number_blocks_in_row(_number_blocks_in_row) {}

	void operator()() const {
		for (int x = 0; x < number_blocks_in_row; x++)
			for (int i = row_start_idx * block_size; i < (row_start_idx + 1)*block_size; i++)
				for (int j = col_start_idx * block_size; j < (col_start_idx + 1)*block_size; j++)
					for (int k = x * block_size; k < (x + 1)*block_size; k++) {
						C[idx(i, j, size)] += A[idx(i, k, size)] * B[idx(k, j, size)];
					}
	}
};

int cannon_tbb(double *&A, double *&B, double *&C, int size, int number_blocks_in_row, task_group &taskGroup) {
	int block_size = size / number_blocks_in_row;
	for (int i = 0; i < number_blocks_in_row; i++)
		for (int j = 0; j < number_blocks_in_row; j++)
			taskGroup.run(CannonTask(A, B, C, i, j, size, block_size, number_blocks_in_row));
	taskGroup.wait();
	return 0;
}

int main(int argc, char** argv) {
	double *A = nullptr;
	double *B = nullptr;
	double *C1 = nullptr;
	double *C2 = nullptr;
	double *C3 = nullptr;
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
	C3 = new double[size*size];
	for (int i = 0; i < size * size; i++) {
		C1[i] = 0;
		C2[i] = 0;
		C3[i] = 0;
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

	const auto start_time_open_mp = omp_get_wtime();
	cannon_par(A, B, C2, size, num_blocks_in_row);
	const auto end_time_open_mp = omp_get_wtime();
	cout << "Cannon OpenMP time    : " << end_time_open_mp - start_time_open_mp << endl;

	task_scheduler_init init(num_thrd);
	task_group task_group;
	auto start_time_tbb = omp_get_wtime();
	cannon_tbb(A, B, C3, size, num_blocks_in_row, task_group);
	auto end_time_tbb = omp_get_wtime();
	cout << "Cannon TBB time       : " << end_time_tbb - start_time_tbb << endl;

	delete[]A;
	delete[]B;
	delete[]C1;
	delete[]C2;
	delete[]C3;

	return 0;
}