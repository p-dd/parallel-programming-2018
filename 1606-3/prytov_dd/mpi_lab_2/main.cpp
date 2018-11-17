#include <mpi.h>
#include <string>
#include <random>
#include <chrono>
#include <iostream>
#include <cmath>
#include <iterator>
#include <vector>
#include <map>

std::vector<int>::iterator binary_search(std::vector<int>& container, int element)
{
	const std::vector<int>::iterator endIt = container.end();

	auto left = container.begin();
	auto right = endIt;

	if (container.empty()
		|| container.front() > element
		|| container.back() < element) {
		return endIt;
	}

	while (distance(left, right) > 0) {
		auto mid = left + distance(left, right) / 2;

		if (element <= *mid) {
			right = mid;
		}
		else {
			left = mid + 1;
		}
	}

	if (*right == element) {
		return right;
	}

	return endIt;
}

std::vector<int>::iterator circular_next(std::vector<int> &l, std::vector<int>::iterator &it)
{
	return std::next(it) == l.end() ? l.begin() : std::next(it);
}
std::vector<int>::iterator circular_prev(std::vector<int> &l, std::vector<int>::iterator &it)
{
	return std::prev(it) == l.begin() - 1 ? l.end() : std::prev(it);
}

bool is_proc_send(int rank, int level, int number, int h_tree)
{
	if (level == 0 || level >= h_tree || rank >= number)
		return false;

	int range = static_cast<int>(pow(2, level)); //4

	if (rank < range)
	{
		if (level == h_tree - 1 && number - range < range)
		{
			return rank < number - range;
		}
		return true;
	}

	return false;
}

bool is_proc_receive(int rank, int level, int number, int h_tree)
{
	if (level == 0 && level >= h_tree)
		return false;
	int range = static_cast<int>(pow(2, level));

	return rank >= range && rank < 2 * range;

}

int MPI_MyScatter(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm communicator)
{
	int proc_rank, proc_num;

	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	MPI_Aint send_extent, recv_extent, send_lb, recv_lb;
	MPI_Type_get_extent(send_datatype, &send_lb, &send_extent);
	MPI_Type_get_extent(recv_datatype, &recv_lb, &recv_extent);

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
		/*if (root == proc_num - 1)
		{
			proc_map.at(root) = 0;
		}*/
		root = 0;
	}
	else
	{
		for (int i = 0; i < proc_num; ++i)
		{
			proc_map.insert({ i, i });	
		}
	}

	std::cout << "Im proc " << proc_map.find(proc_rank)->first << " of " << proc_num << " my old rank " << proc_map.find(proc_rank)->second << std::endl;
	//dedug
	MPI_Barrier(communicator);

	double h = log2(proc_num);
	int height = static_cast<int>(h);
	if (h - height > 0)
		++height;

	//debug
	if (proc_rank == root)
		std::cout <<proc_rank << " tree height: " << height << std::endl << std::endl;
	
	char* buffer = nullptr;
	int buffer_start = 0;

	if (proc_rank == root)
	{
		buffer = static_cast<char*>(send_data);
		if (buffer == nullptr)
			std::cout << "ERROR!" << std::endl;
	}

	std::vector<int> circular_proc_set(proc_num);

	circular_proc_set.reserve(proc_num);
	for (int i = 0; i < proc_num; ++i)
		circular_proc_set.push_back(i);

	auto destination_proc = circular_proc_set.begin();

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
		int level_buffer_size = 0;
		// блок выделения памяти на каждой итерации
		// проверить что будет происходить когда число процессов будет равно 5? 7? 9? 11?

		int size = (send_count * proc_num * send_extent) / (2 * (level_tree + 1));
		if (size < send_count * send_extent)
			size = send_count * send_extent;

		level_buffer = new char[size];
		if (level_tree == 0)
		{
			int tmp = proc_num / 2;
			//if (tmp % 2 == 1)
			//	--tmp;
			level_buffer_size = send_count * tmp;
		}

		if (level_tree > 0)
		{
			//level_buffer_size = proc_num * send_count / ((level_tree + 1) * 2); // 4*2/4 = 4 -> its problem on level 0 when proc_num == 3
			//if (level_tree == 2)
			//{
			//	std::cout << "_________________________hello from " << proc_rank << " at lvl " << level_tree << std::endl;
			//	MPI_Barrier(communicator);
			//}
			if (is_proc_send(proc_rank + 2, level_tree + 1, proc_num, height) || (is_proc_receive(proc_rank, level_tree, proc_num, height) && is_proc_send(proc_rank, level_tree + 1, proc_num, height)))
			{
				level_buffer_size = send_count * (height - level_tree);
				if (level_tree == 2 && proc_rank == 2)
					std::cout << "_________________________buffer != send count" << std::endl;
				//std::cout << ">>>>>>>>>>lvl " << level_tree << ": " << proc_rank << ", is proc send: " << is_proc_send(proc_num, level_tree, proc_num, height) << ", cut size " << " level buffer size: " << level_buffer_size << std::endl << std::endl;
			}
			else
			{
				if (level_tree == 2 && proc_rank == 2)
					std::cout << "_________________________buffer = send count" << std::endl;
				level_buffer_size = send_count;
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		//пока так в дальнейшем может быть из-за этого ошибка, обратить внимание !!!!!//
		////////////////////////////////////////////////////////////////////////////////

		if (level_buffer_size < send_count)
			level_buffer_size = send_count;

		//}
		// блок посылки первого сообщения от root к root + 1

		if (proc_rank == 0 && level_tree == 0) // go at level 0
		{
			std::cout << "lvl " << level_tree << ": " << proc_rank << " root send data to " << proc_map.find(1)->first << " level buffer size: " << level_buffer_size << std::endl << std::endl;

			if (buffer == nullptr)
				std::cout << "ERROR! 1" << std::endl;

			MPI_Send(buffer, level_buffer_size, send_datatype, proc_map.at(1), 0, communicator);

			buffer_start = level_buffer_size;

			for (int i = (proc_num - 1) * recv_extent * recv_count, j = 0; i < proc_num * recv_count * recv_extent; ++i, ++j)
			{
				*(static_cast<char*>(recv_data) + j) = *(static_cast<char*>(send_data) + i);
			}
			continue;
		}

		if (proc_rank == 1 && level_tree == 0) // go at level 0
		{
			//std::cout << "lvl " << level_tree << ": " << proc_rank << " msg recv from root: " << proc_map.at(0) << std::endl << std::endl;

			if (level_buffer == nullptr)
				std::cout << "ERROR! 2" << std::endl;

			MPI_Recv(level_buffer, level_buffer_size, recv_datatype, proc_map.at(0), 0, communicator, MPI_STATUS_IGNORE);

			std::cout << "lvl " << level_tree << ": " << proc_rank << " msg recv from root: " << root << std::endl << std::endl;

			if (!is_proc_send(proc_rank, level_tree + 1, proc_num, height)/*level_tree + 1 == height || (proc_num % 2 == 1 && proc_num == 3)*/)
			{
				for (int i = 0, j = 0; i < level_buffer_size * recv_extent; ++i, ++j)
				{
					*(static_cast<char*>(recv_data) + j) = *(level_buffer + i);
				}
				continue;
			}
			if (level_buffer_size % (send_count * 2) == 0)
			{
				for (int i = recv_extent * level_buffer_size / 2, j = 0; i < level_buffer_size * recv_extent; ++i, ++j)
				{
					*(static_cast<char*>(recv_data) + j) = *(level_buffer + i);
				}
				std::cout << "``1``````````lvl " << level_tree << ": " << proc_rank << " HELLO I Write MY DATA" << std::endl << std::endl;
			}
			else
			{
				for (int i = recv_extent * (level_buffer_size - send_count), j = 0; i < level_buffer_size * recv_extent; ++i, ++j)
				{
					*(static_cast<char*>(recv_data) + j) = *(level_buffer + i);
				}
				//for (int i = 0; i < send_count * recv_extent; ++i)
				//{
				//	*(static_cast<char*>(recv_data) + i) = *(level_buffer + i);
				//}
				std::cout << "``2``````````lvl " << level_tree << ": " << proc_rank << " HELLO I Write MY DATA" << std::endl << std::endl;
			}
			buffer = level_buffer;


			/*for (int i = recv_extent * level_buffer_size / 2, j = 0; i < level_buffer_size * recv_extent; ++i, ++j)
			{
				*(static_cast<char*>(recv_data) + j) = *(level_buffer + i);
			}
			buffer = level_buffer;*/
			std::cout << "lvl " << level_tree << ": " << proc_rank << " write buffer" << std::endl << std::endl;
			
			continue;
		}

		// блок дальнейшей рассылки сообщений между процессами с периодом 

		if (is_proc_send(proc_rank, level_tree, proc_num, height))
			///*proc_rank == root || */proc_rank <= static_cast<int>(pow(2, level_tree)) /*root + level_tree*/ && proc_num % 2 == 0
			//|| proc_rank/*proc_rank == root + level_tree &&*/ /*proc_rank != proc_num - 1 &&*/ proc_num % 2 == 1 /*&& proc_num - 1 != proc_rank + level_tree*/)
												 // когда больше 4х процессов написать условие сложнее: чтобы входили все процессы ранг который лежит от рута до root + level tree
		{													// третье условие говорит срабатавает когда нечетное количество процессов и нужно отправить в конце только с рута
			
			if (proc_rank == root)
			{
				//std::cout << "lvl " << level_tree << ": " << proc_rank << " root send data to " << proc_map.find(proc_rank + level_tree * 2)->first << " level buffer size: " << level_buffer_size << std::endl << std::endl;
				MPI_Send(buffer + buffer_start * send_extent, level_buffer_size, send_datatype, proc_map.at(proc_rank + level_tree * 2), 0, communicator);

				buffer_start += level_buffer_size;

				std::cout << "lvl "<< level_tree <<": " << proc_rank << " root send data to " << proc_map.find(proc_rank + level_tree * 2)->first << " level buffer size: " << level_buffer_size << std::endl << std::endl;
				continue;
			}
			//std::cout << "lvl " << level_tree << ": " << proc_rank << " ne root send data to " << proc_map.find(proc_rank + level_tree * 2)->first << " level buffer size: " << level_buffer_size << std::endl << std::endl;
			MPI_Send(buffer, level_buffer_size, send_datatype, proc_map.at(proc_rank + level_tree * 2), 0, communicator);
			std::cout << "lvl " << level_tree << ": " << proc_rank << " ne root send data to "<< proc_map.find(proc_rank + level_tree * 2)->first << " level buffer size: " << level_buffer_size << std::endl << std::endl;
			continue;

		}

		if (is_proc_receive(proc_rank, level_tree, proc_num, height))
		{
			//std::cout << "------------------lvl " << level_tree << ": " << proc_rank << " msg receive from " << proc_rank - 2 * level_tree << std::endl << std::endl;
			MPI_Recv(level_buffer, level_buffer_size, recv_datatype, proc_map.at(proc_rank - 2 * level_tree), 0, communicator, MPI_STATUS_IGNORE);

			std::cout << "------------------lvl " << level_tree << ": " << proc_rank << " msg receive from " << proc_rank - 2 * level_tree << " test info: "<< is_proc_send(proc_rank, level_tree + 1, proc_num, height) << std::endl << std::endl;

			if (!is_proc_send(proc_rank, level_tree + 1, proc_num, height)/*(proc_num % 2 == 1 && proc_rank > proc_num / 2 - 1) 
				|| (proc_num % 2 == 0 && proc_rank > (proc_num / 2) - 1) && (proc_num / 2) % 2 == 0 || (proc_num % 2 == 0 && proc_rank > (proc_num / 2) - 2) && (proc_num / 2) % 2 == 1*/)
			{
				for (int i = 0; i < level_buffer_size * recv_extent; ++i)
				{
					*(static_cast<char*>(recv_data) + i) = *(level_buffer + i);
				}
				std::cout << "lvl " << level_tree << ": " << proc_rank << " bye :3" << std::endl << std::endl;
			}
			else
			{
				if (level_buffer_size % 2 == 0)
				{
					for (int i = recv_extent * level_buffer_size / 2, j = 0; i < level_buffer_size * recv_extent; ++i, ++j)
					{
						*(static_cast<char*>(recv_data) + j) = *(level_buffer + i);
					}
					std::cout << "````be here````````lvl " << level_tree << ": " << proc_rank << " HELLO I Write MY DATA" << std::endl;
				}
				else
				{
					for (int i = recv_extent * (level_buffer_size - send_count), j = 0; i < level_buffer_size * recv_extent; ++i, ++j)
					{
						*(static_cast<char*>(recv_data) + j) = *(level_buffer + i);
					}
					//for (int i = 0; i < send_count * recv_extent; ++i)
					//{
					//	*(static_cast<char*>(recv_data) + i) = *(level_buffer + i);
					//}
					std::cout << "````````````lvl " << level_tree << ": " << proc_rank << " HELLO I Write MY DATA" << std::endl << std::endl;
				}
				if (level_tree == 1)
				{
					std::cout << "------>>>>>hello from " << proc_rank << " at lvl " << level_tree << std::endl;
				}
				buffer = level_buffer;
				if (level_tree == 1)
				{
					std::cout << "------>>>>>hello from " << proc_rank << " at lvl " << level_tree << std::endl;
				}
			}
			if (level_tree == 1)
			{
				std::cout << "-defore delete----->>>>>hello from " << proc_rank << " at lvl " << level_tree << std::endl;
			}
			//delete[] level_buffer;
			if (level_tree == 1)
			{
				std::cout << "-after delete----->>>>>hello from " << proc_rank << " at lvl " << level_tree << std::endl;
			}
		}
	}

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