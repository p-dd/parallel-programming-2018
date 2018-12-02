#include <mpi.h>
#include <string>
#include <random>
#include <chrono>
#include <iostream>
#include <iterator>
#include <vector>
#include <map>

int MPI_MyScatter(void* send_buf, int send_count, MPI_Datatype send_type, void* recv_buf, int recv_count, MPI_Datatype recv_type, int root, MPI_Comm comm)
{
	int proc_rank, proc_num;
	MPI_Comm_size(comm, &proc_num);
	int error = MPI_SUCCESS;
	MPI_Aint send_extent;
	MPI_Aint recv_extent;
	MPI_Comm_rank(comm, &proc_rank);
	MPI_Aint send_lb, recv_lb;
	MPI_Type_get_extent(send_type, &send_lb, &send_extent);
	MPI_Type_get_extent(recv_type, &recv_lb, &recv_extent);
	MPI_Status status;

	char * tree_buf;

	std::map<int, int> proc_map;
	
	if (root != 0)
	{
		int shift = root;
		if (proc_rank < root)
			proc_rank = proc_rank - shift + proc_num;
		else
			proc_rank -= shift;
		for (int i = 0; i < proc_num; ++i)
		{
			if (i <= proc_num - 1 - root)
				proc_map.insert({ i, i + shift });
			else
				proc_map.insert({ i, i - proc_num + shift }); // >root
		}
		root = 0;
	}
	else 
	{
		for (int i = 0; i < proc_num; ++i)
		{
			proc_map.insert({ i, i });	
		}
	}
	//std::cout << "Im proc " << proc_map.find(proc_rank)->first << " of " << proc_num << " my old rank " << proc_map.find(proc_rank)->second << std::endl;
	//MPI_Barrier(comm);

	if (proc_map.at(proc_rank) % 2)
	{
		error = MPI_Recv(recv_buf, recv_count, recv_type, proc_map.at(proc_rank - 1), 0, comm, &status);
	}
	else
	{
		int tree_part = proc_num;
		int x = proc_map.at(proc_rank);
		while (x % tree_part != 0)
		{
			x %= tree_part;
			tree_part >>= 1;
		}
		if (proc_map.at(proc_rank) != root)
		{
			tree_buf = new char[tree_part * recv_extent * recv_count];
			error = MPI_Recv(tree_buf, tree_part * recv_count, recv_type, proc_map.at(proc_rank - tree_part), 0, comm, &status);
		}
		else
		{
			tree_buf = static_cast<char*>(send_buf);
		}

		for (int i = 0; i < send_count * send_extent; i++) 
		{
			*(static_cast<char*>(recv_buf) + i) = *(tree_buf + i);
		}

		tree_part >>= 1;

		for (tree_part; tree_part >= 1; tree_part >>= 1) 
		{
			error = MPI_Send(tree_buf + tree_part * send_count * send_extent, tree_part * send_count, send_type, proc_map.at(tree_part + proc_rank), 0, comm);
		}
	}
	return error;
}

int main(int argc, char* argv[])
{
	int rank, proc;
	double* sys_mas = nullptr;
	double* my_mas = nullptr;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc);

	const int size = stoi(std::string(argv[1]));

	const int root = stoi(std::string(argv[2]));

	if (rank == root)
	{
		std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());

		const std::uniform_real_distribution <double> distribution(-100, 100);

		sys_mas = new double[size];
		my_mas = new double[size];
		if (size < 21) {
			std::cout.precision(4);
			for (int i = 0; i < size; ++i)
			{
				sys_mas[i] = distribution(generator);
				my_mas[i] = sys_mas[i];
			}
			std::cout << std::endl;
			std::cout << "sys_mas : ";
			for (int i = 0; i < size; ++i)
			{
				std::cout << sys_mas[i] << " ";
			}
			std::cout << std::endl << "my_mas : ";
			for (int i = 0; i < size; ++i)
			{
				std::cout << my_mas[i] << " ";
			}
			std::cout << std::endl;
			std::cout << std::endl;
		}
	}


	const int part = size / proc;
	double* sys_arr = new double[part];
	double* my_arr = new double[part];
	double start_time, end_time;
	if (rank == root)
	{
		start_time = MPI_Wtime();
	}

	MPI_Scatter(sys_mas, part, MPI_DOUBLE, sys_arr, part, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == root)
	{
		end_time = MPI_Wtime();
		std::cout << " MPI_Scatter time: " << end_time - start_time << std::endl;
	}
	
	if (rank == root)
	{
		start_time = MPI_Wtime();
	}

	MPI_MyScatter(my_mas, part, MPI_DOUBLE, my_arr, part, MPI_DOUBLE, root, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == root)
	{
		end_time = MPI_Wtime();
		std::cout << " MPI_MyScatter time: " << end_time - start_time << std::endl << std::endl;
	}
	
	/*std::cout.precision(4);
	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "I`m " << rank << " proc my sys_data: ";
	for (int i = 0; i < part; ++i)
		std::cout << sys_arr[i] << " ";*/

	MPI_Barrier(MPI_COMM_WORLD);
	std::cout.precision(4);
	if (size < 21) {
		std::cout << "I`m " << rank << " proc my my_data: ";
		for (int i = 0; i < part; ++i)
			std::cout << my_arr[i] << " ";
		std::cout << std::endl;
	}
	MPI_Finalize();

	return 0;
}