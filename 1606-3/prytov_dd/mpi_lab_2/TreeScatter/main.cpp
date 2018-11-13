#include <mpi.h>
#include <string>
#include <random>
#include <chrono>
#include <iostream>
#include <cmath>
#include <iterator>

int calculate_source()
{
	return 0;
}

int calculate_destination()
{
	return 0;
}


std::vector<int>::iterator circular_next(std::vector<int> &l, std::vector<int>::iterator &it)
{
	return std::next(it) == l.end() ? l.begin() : std::next(it);
}
std::vector<int>::iterator circular_prev(std::vector<int> &l, std::vector<int>::iterator &it)
{
	return std::prev(it) == l.begin() - 1 ? l.end() : std::prev(it);
}


/*struct node
{
	int process_num = 0;
	int destination_proc_num = 0;
};

class cyclic_list
{
public:
	enum status { receive, send };
	cyclic_list *next;

	cyclic_list(int num = 0, cyclic_list* next = nullptr) : process_num(num), next(next) {}
	void set_next(cyclic_list* cl) 
	{
		next = cl;
	}
	cyclic_list* get_next() const
	{
		return next;
	}
};*/

int MPI_MyScatter(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm communicator)
{
	int proc_rank, proc_num;

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	std::cout << "Im proc " << proc_rank << " of " << proc_num << std::endl;

	MPI_Barrier(communicator);
	MPI_Aint send_extent, recv_extent, send_lb, recv_lb;
	MPI_Type_get_extent(send_datatype, &send_lb, &send_extent);
	MPI_Type_get_extent(recv_datatype, &recv_lb, &recv_extent);
	                                              //MPI_Recv(buffer, count, datatype, _from(_r / 2, root, procNum), 0, comm, MPI_STATUS_IGNORE); use it!
	int height = log(proc_num) + 1;

	//debug
	if (proc_rank == root)
		std::cout << "tree height: " << height << std::endl;
	
	char* buffer = nullptr;
	int buffer_start = 0;

	if (proc_rank == root)
	{
		buffer = static_cast<char*>(send_data);
	}

	/*std::vector<int> circular_proc_set(proc_num);

	circular_proc_set.reserve(proc_num);
	for (int i = 0; i < proc_num; ++i)
		circular_proc_set.push_back(i);

	auto destination_proc = circular_proc_set.begin();*/

	if (proc_num == 1)
	{
		for (int i = 0; i < recv_count * recv_extent; ++i)
		{
			*(static_cast<char*>(recv_data) + i) = *(static_cast<char*>(send_data) + i);
		}
		return 0;
	}

	for (int level_tree = 0; level_tree < height; ++level_tree) 
	{
		char* level_buffer;
		int level_buffer_size;

		if((proc_num % 2 == 1) && (level_tree == 0))
		{
			level_buffer = new char[((send_count * proc_num * send_extent) / 2) - (send_extent / 2)]; 

			level_buffer_size = proc_num * send_count / 2; 

		}
		else
		{
			int tmp = proc_num;
			if (proc_num % 2 == 1)
				++tmp;
		
			level_buffer = new char[(send_count * tmp * send_extent) / (2 * (level_tree + 1))]; //4*8 / 4 = 8 +

			level_buffer_size = tmp * send_count / ((level_tree + 1) * 2); // 4*1 / 4 = 1 +
		}

		if (proc_rank == root && level_tree == 0) // go at level 0
		{
			if (proc_num % 2 == 1)
			{
				MPI_Send(buffer, level_buffer_size, send_datatype, (proc_rank + 1) % proc_num, 0, communicator);

				buffer_start = level_buffer_size;
			}
			else
			{
				MPI_Send(buffer, level_buffer_size, send_datatype, (proc_rank + 1) % proc_num, 0, communicator);

				buffer_start = level_buffer_size;
			}
			for (int i = (proc_num - 1) * recv_extent * recv_count, j = 0; i < proc_num * recv_count * recv_extent; ++i, ++j)
			{
				*(static_cast<char*>(recv_data) + j) = *(static_cast<char*>(send_data) + i);
			}
			continue;
		}

		if (proc_rank == root || ((proc_rank == ((root + level_tree) % proc_num)) && (proc_rank % 2 == 0))
			|| ((proc_rank == ((root + level_tree) % proc_num)) && (root + level_tree * 2 != proc_num - 1) && (proc_rank % 2 == 1))) // когда больше 4х процессов написать условие сложнее: чтобы входили все процессы ранг который лежит от рута до root + level tree
		{																													     // третье условие говорит срабатавает когда нечетное количество процессов и нужно отправить в конце только с рута
			
			if (proc_rank == root)
			{ //необходимо сделать сдвижение на дельту которая накапливается с течением времени
				MPI_Send(buffer + buffer_start * send_extent, level_buffer_size, send_datatype, (proc_rank + level_tree * 2) % proc_num, 0, communicator);

				buffer_start += level_buffer_size;

				std::cout << "lvl 1: " << proc_rank << " root send data to " << (proc_rank + level_tree * 2) % proc_num << std::endl;
				continue;
			}

			MPI_Send(buffer, level_buffer_size, send_datatype, (proc_rank + level_tree * 2) % proc_num, 0, communicator);
			std::cout << "lvl 1: " << proc_rank << " ne root send data to "<< (proc_rank + level_tree * 2) % proc_num << std::endl;
			continue;

		}

		if ((proc_rank == ((root + 1) % proc_num)) && (level_tree == 0)) // go at level 0
		{
			MPI_Recv(level_buffer, level_buffer_size, recv_datatype, root, 0, communicator, MPI_STATUS_IGNORE);

			std::cout << "lvl 0: " << proc_rank << " msg recv from root: " << root << std::endl;

			if (level_tree + 1 == height || (proc_num % 2 == 1))
			{
				for (int i = 0, j = 0; i < level_buffer_size * recv_extent; ++i, ++j)
				{
					*(static_cast<char*>(recv_data) + j) = *(level_buffer + i);
				}
			}
			else
			{
				for (int i = recv_extent * level_buffer_size / 2, j = 0; i < level_buffer_size * recv_extent; ++i, ++j)
				{
					*(static_cast<char*>(recv_data) + j) = *(level_buffer + i);
				}
				buffer = level_buffer;
				std::cout << "lvl 0: " << proc_rank << " write buffer" << std::endl;
			}
			continue;
		}

		if (proc_rank == ((root + level_tree * 2) % proc_num) || ((proc_rank == (root + level_tree + level_tree * 2) % proc_num) && (proc_num % 2 == 0)) ) // частный случай когда 4 процесса 
		{

			int destination = proc_rank - 2;

			if(destination < 0)
			{
				if (proc_rank % 2 == 1)
					destination = proc_num - 1;
				else
					destination = proc_num - 2;
			}

			MPI_Recv(level_buffer, level_buffer_size, recv_datatype, destination, 0, communicator, MPI_STATUS_IGNORE);

			std::cout << "------------------lvl 1: " << proc_rank << " msg receive from " << destination << std::endl;
			/*auto test_buf = reinterpret_cast<double*>(level_buffer);
			if (proc_rank == 2)
			{
				std::cout << "data test: " << std::endl;
				for (int i = 0; i < level_buffer_size; ++i)
				{
					std::cout << *(test_buf + i) << std::endl;
				}
			}*/
			//std::cout << "lvl 1: " << proc_rank << " msg receive from " << proc_rank - 2 << std::endl;

			if (level_tree + 1 == height)
			{
				for (int i = 0; i < level_buffer_size * recv_extent; ++i)
				{
					*(static_cast<char*>(recv_data) + i) = *(level_buffer + i);
				}
				std::cout << "lvl 1: " << proc_rank << " bye :3" << std::endl;
			}
			else
			{
				for (int i = recv_extent * level_buffer_size / 2, j = 0; i < level_buffer_size * recv_extent; ++i, ++j)
				{
					*(static_cast<char*>(recv_data) + j) = *(level_buffer + i);
				}
				buffer = level_buffer;
			}
		}
	}

	/*if (buffer != nullptr)
	{
		for (int i = 0; i < recv_count * recv_extent; ++i)
		{
			*(static_cast<char*>(recv_data) + i) = *(buffer + i);
		}
	}*/

	/*if (proc_rank == root)
	{
		for (int i = (proc_num - 1) * recv_extent * recv_count, j = 0; i < proc_num * recv_count * recv_extent; ++i, ++j)
		{
			*(static_cast<char*>(recv_data) + j) = *(static_cast<char*>(send_data) + i);
		}
	}*/
	return 0;
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


	const int part = size / proc;
	double* sys_arr = new double[part];
	double* my_arr = new double[part];

	MPI_Scatter(sys_mas, part, MPI_DOUBLE, sys_arr, part, MPI_DOUBLE, root, MPI_COMM_WORLD);

	MPI_MyScatter(my_mas, part, MPI_DOUBLE, my_arr, part, MPI_DOUBLE, root, MPI_COMM_WORLD);

	std::cout.precision(4);
	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "I`m " << rank << " proc my sys_data: ";
	for (int i = 0; i < part; ++i)
		std::cout << sys_arr[i] << " ";

	/*std::cout << "global data in " << rank << " proc: ";
	for (int i = 0; i < size; ++i)
		std::cout << sys_mas[i] << " ";
	std::cout << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);*/

	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "I`m " << rank << " proc my my_data: ";
	for (int i = 0; i < part; ++i)
		std::cout << my_arr[i] << " ";

	std::cout << std::endl;

	MPI_Finalize();

	return 0;
}