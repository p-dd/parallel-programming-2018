#include <mpi.h>
#include <cstdlib>
#include <time.h>
#include <Windows.h>
#include <queue>
#include <iostream>
#include <string>

#define INFINITI 10000000
#define WEIGHT 5
using namespace std;

void print_d(int* d, int size) {
	for (int i = 0; i < size; i++) {
		cout << d[i] << " ";
	}
	cout << endl;
}
void print_graph(int** G, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << G[i][j] << " ";
		}
		cout << endl;
	}
}

int* init_graph(int countEdge, int countVertex) {
	int rank, procNum;

	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int* one_graph = new int[countVertex * countVertex];
	int** graph = new int*[countVertex];

	for (int i = 0; i < countVertex; i++) {
		graph[i] = new int[countVertex];
	}
	srand(time(NULL));

	for (int i = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			if (i == j) {
				graph[i][j] = 0;

			}
			else {
				graph[i][j] = rand() % 3;
				if (graph[i][j] == 0) {
					graph[i][j] = INFINITI;
				}
			}
			//cout << G[i][j] << " ";
		}
		//cout << endl;
	}
	for (int i = 0, t = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			one_graph[t] = graph[j][i];
			t++;
		}
	}
	return one_graph;
}

int* init_d(int size) {
	int* d = new int[size];
	for (int i = 0; i < size; i++) {
		d[i] = INFINITI;
	}
	return d;
}

int* dijkstra(int* graph, int start, int count_vertex) {
	int* d = init_d(count_vertex);
	d[start] = 0;
	priority_queue<pair<int, int>>  queue;
	queue.push(make_pair(0, start));
	while (!queue.empty()) {
		int v = queue.top().second, cur_d = queue.top().first;
		queue.pop();

		if (cur_d > d[v])  continue;

		for (int i = 0; i < count_vertex; ++i) {
			int to = i,	len = graph[i * count_vertex + v];
			if (d[v] + len < d[to]) {
				d[to] = d[v] + len;
				queue.push(make_pair(d[to], to));
			}
		}
	}
	return d;
}

int* parallel_dijkstra(int* graph, int start, int count_vertex) {
	int* d = nullptr;
	
	int rank = 0, proc_num = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int* dist_d = new int[proc_num];
	int* dist_graph = new int[proc_num];
	int* count_element_d = new int[proc_num];
	int* count_element_graph = new int[proc_num];
	int part_size = count_vertex / (proc_num - 1);
	
	count_element_d[0] = 0; //на нулевой не шлем 
	dist_d[0] = 0;
	count_element_graph[0] = 0;
	dist_graph[0] = 0;
	for (int i = 1; i < proc_num - 1; i++) {
		count_element_d[i] = part_size;
		dist_d[i] = (i - 1) * part_size;
		count_element_graph[i] = count_vertex * part_size;
		dist_graph[i] = (i - 1) * count_vertex * part_size;
	}
	dist_d[proc_num - 1] = (proc_num - 2) * part_size;
	count_element_d[proc_num - 1] = count_vertex - (proc_num - 2) * part_size;
	dist_graph[proc_num - 1] = (proc_num - 2) * (count_vertex * part_size);
	count_element_graph[proc_num - 1] = count_vertex * count_vertex - (proc_num - 2) * (count_vertex*part_size);

	int* part_graph = new int[count_element_graph[rank]];
	int* part_d = new int[count_element_d[rank]];
	int flag = 1, flag_finish_find_current_min_distation_in_sigment_of_rank = 1;

	pair<int, int> current_vertex;
	pair<int, int> current_vertex_with_current_min_destation;

	MPI_Status st;
	if (rank == 0) {
		MPI_Scatterv(graph, count_element_graph, dist_graph, MPI_INT, part_graph, count_element_graph[rank], MPI_INT, 0, MPI_COMM_WORLD); //пересылка матрицы смежности
		d = init_d(count_vertex);
		d[start] = 0;
		priority_queue<pair<int, int>>  queue;
		queue.push(make_pair(0, start));

		int flag_rank_0 = 1;
		while (!queue.empty()) {

			int v = queue.top().second;
			int	dest = queue.top().first;
			queue.pop();
			if (dest > d[v])
				continue;
			MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
			current_vertex = make_pair(dest, v);
			MPI_Bcast(&current_vertex, 1, MPI_2INT, 0, MPI_COMM_WORLD); //рассылка выбранной вершины и массива длин путей
			MPI_Scatterv(d, count_element_d, dist_d, MPI_INT, part_d, count_element_d[rank], MPI_INT, 0, MPI_COMM_WORLD); 
			for (int i = 1; i < proc_num; i++) {
				flag_rank_0 = 1;
				while (flag_rank_0 != 0) {

					MPI_Recv(&flag_finish_find_current_min_distation_in_sigment_of_rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &st);
					flag_rank_0 = flag_finish_find_current_min_distation_in_sigment_of_rank;
					if (flag_rank_0 == 0) {
						continue;
					}
					MPI_Recv(&current_vertex_with_current_min_destation, 1, MPI_2INT, i, 0, MPI_COMM_WORLD, &st);
					d[current_vertex_with_current_min_destation.second] = current_vertex_with_current_min_destation.first;

					queue.push(make_pair(current_vertex_with_current_min_destation.first, current_vertex_with_current_min_destation.second));
				}
			}
			flag_rank_0 = 1;
		}
		flag = 0;
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else {
		MPI_Scatterv(graph, count_element_graph, dist_graph, MPI_INT, part_graph, count_element_graph[rank], MPI_INT, 0, MPI_COMM_WORLD); //прием блока матрицы смежности

		while (flag) {
			MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (flag == 0) {
				continue;
			}
			MPI_Bcast(&current_vertex, 1, MPI_2INT, 0, MPI_COMM_WORLD); //прием стартовой вершины и массива длин путей
			MPI_Scatterv(d, count_element_d, dist_d, MPI_INT, part_d, count_element_d[rank], MPI_INT, 0, MPI_COMM_WORLD);

			for (int j = 0; j < count_element_d[rank]; j++) {

				int to = j, len = part_graph[j * count_vertex + current_vertex.second];

				if (current_vertex.first + len < part_d[to]) {
					flag_finish_find_current_min_distation_in_sigment_of_rank = 1;
					MPI_Send(&flag_finish_find_current_min_distation_in_sigment_of_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
					current_vertex_with_current_min_destation.first = current_vertex.first + len;
					current_vertex_with_current_min_destation.second = to + dist_d[rank];

					MPI_Send(&current_vertex_with_current_min_destation, 1, MPI_2INT, 0, 0, MPI_COMM_WORLD);
				}
			}

			flag_finish_find_current_min_distation_in_sigment_of_rank = 0;
			MPI_Send(&flag_finish_find_current_min_distation_in_sigment_of_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
	return d;
}

bool is_correct_implementation(int* d1, int* d2, int size) {
	int count_mistakes = 0;
	for (int i = 0; i < size; i++) {
		if (d1[i] != d2[i]) {
			count_mistakes++;
			return false;
		}
	}
	if (count_mistakes == 0) {
		return true;
	}
}
int* init_route(int countVertex) {
	int** graph = new int*[countVertex];
	int* one_graph = new int[countVertex*countVertex];
	for (int i = 0; i < countVertex; i++) {
		graph[i] = new int[countVertex];
	}
	srand(time(NULL));
	for (int i = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			if (i == j) {
				graph[i][j] = 0;
			}
			else {
				if (j == i + 1) {
					graph[i][j] = 1 + rand() % WEIGHT;
				}
				else {
					graph[i][j] = INFINITI;
				}
			}
		}
	}
	for (int i = 0, t = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			one_graph[t] = graph[j][i];
			t++;
		}
	}
	return one_graph;
}

int* ideal_d_for_route(int* graph, int start, int count_vertex) {
	int* d = init_d(count_vertex);
	d[start] = 0;
	for (int ofset = 0, i = 1 + start; (i < count_vertex); ofset++, i++) {
		d[i] = d[i - 1] + graph[i*count_vertex + start + ofset];
	}
	return d;
}

bool test_route(int count_vertex) {
	int start = rand() % count_vertex;
	int* graph = init_route(count_vertex);
	int* d = dijkstra(graph, start, count_vertex);
	int* d_ideal = ideal_d_for_route(graph, start, count_vertex);
	return is_correct_implementation(d, d_ideal, count_vertex);
}

int* init_round(int count_vertex) {
	int** graph = new int*[count_vertex];
	int* one_graph = new int[count_vertex*count_vertex];
	for (int i = 0; i < count_vertex; i++) {
		graph[i] = new int[count_vertex];
	}
	srand(time(NULL));
	for (int i = 0; i < count_vertex; i++) {
		for (int j = 0; j < count_vertex; j++) {
			if (i == j) {
				graph[i][j] = 0;
			}
			else {
				if (j == i + 1) {
					graph[i][j] = 1 + rand() % WEIGHT;
				}
				else {
					graph[i][j] = INFINITI;
				}
			}
		}
	}
	graph[count_vertex - 1][0] = 1 + rand() % WEIGHT;
	for (int i = 0, t = 0; i < count_vertex; i++) {
		for (int j = 0; j < count_vertex; j++) {
			one_graph[t] = graph[j][i];
			t++;
		}
	}
	return one_graph;
}

int* ideal_d_for_round(int* graph, int start, int count_vertex) {
	int* d = init_d(count_vertex);
	d[start] = 0;
	for (int ofset = 0, i = 1 + start; (i < count_vertex); ofset++, i++) {
		d[i] = d[i - 1] + graph[i*count_vertex + start + ofset];
	}
	for (int i = 0; i < start; i++) {
		if (i == 0) {
			d[i] = d[count_vertex - 1] + graph[count_vertex - 1];
		}
		else {
			d[i] = d[i - 1] + graph[i*count_vertex + (i - 1)];
		}
	}
	return d;
}

bool test_round(int count_vertex) {
	int start = rand() % count_vertex;
	int* graph = init_round(count_vertex);
	int* d = dijkstra(graph, start, count_vertex);
	int* dIdeal = ideal_d_for_round(graph, start, count_vertex);
	return is_correct_implementation(d, dIdeal, count_vertex);
}

int* init_star(int count_vertex) {
	int** G = new int*[count_vertex];
	int* oneG = new int[count_vertex*count_vertex];
	for (int i = 0; i < count_vertex; i++) {
		G[i] = new int[count_vertex];
	}
	srand(time(NULL));
	for (int i = 0; i < count_vertex; i++) {
		for (int j = 0; j < count_vertex; j++) {
			if (i == j) {
				G[i][j] = 0;
			}
			else {
				if (i == 0) {
					G[i][j] = 1 + rand() % WEIGHT;
				}
				else {
					G[i][j] = INFINITI;
				}
			}
		}
	}
	for (int i = 0, t = 0; i < count_vertex; i++) {
		for (int j = 0; j < count_vertex; j++) {
			oneG[t] = G[j][i];
			t++;
		}
	}
	return oneG;
}

int* ideal_d_for_star(int* graph, int start, int count_vertex) {
	int* d = init_d(count_vertex);
	d[start] = 0;
	for (int i = 0; i < count_vertex; i++) {
		d[i] = graph[i*count_vertex + start];
	}
	return d;
}

bool test_star(int count_vertex) {
	int start = rand() % count_vertex;
	int* graph = init_star(count_vertex);
	int* d = dijkstra(graph, start, count_vertex);
	int* d_ideal = ideal_d_for_star(graph, start, count_vertex);
	return is_correct_implementation(d, d_ideal, count_vertex);
}

bool test_for_iter_algoritm(int count_vertex) {
	bool test = true;
	cout << "testing step by step implementation: " << endl;
	if (!test_route(count_vertex)) {
		cout << "	test route is falid" << endl;
		test = false;
	}
	else {
		cout << "	test route is sucssed" << endl;
	}
	if (!test_round(count_vertex)) {
		cout << "	test round is falid" << endl;
		test = false;
	}
	else {
		cout << "	test round is sucssed" << endl;
	}
	if (!test_star(count_vertex)) {
		cout << "	test star is falid" << endl;
		test = false;
	}
	else {
		cout << "	test star is sucssed" << endl;
	}
	return test;
}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	const int count_vertex = stoi(std::string(argv[1]));
	srand(time(NULL));
	const int count_edge = (count_vertex - 1) + rand() % ((count_vertex * (count_vertex - 1)) / 2);

	int* d_step = nullptr;
	int* d_parallel = nullptr;
	int* graph = nullptr;

	double t1_step = 0.0, t2_step = 0.0, t1_parallel = 0.0, t2_parallel = 0.0;

	int start = rand() % (count_vertex - 1);

	if (proc_rank == 0) {
		graph = init_graph(count_edge, count_vertex);
	}

	if (proc_rank == 0) {
		t1_parallel = MPI_Wtime();
	}
	d_parallel = parallel_dijkstra(graph, start, count_vertex);
	MPI_Barrier(MPI_COMM_WORLD);
	if (proc_rank == 0) {
		t2_parallel = MPI_Wtime();

		cout << endl << " Parallel Dijkstra time: " << t2_parallel - t1_parallel << endl << endl;

		/*for (int i = 0; i < count_vertex; i++) {
			cout << d_parallel[i] << " ";
		}
		cout << endl << endl;*/
		t1_step = MPI_Wtime();
		d_step = dijkstra(graph, start, count_vertex);
		t2_step = MPI_Wtime();
		cout << endl << " Dijkstra time: " << t2_step - t1_step << endl << endl;
		/*for (int i = 0; i < count_vertex; i++) {
			cout << d_step[i] << " ";
		}
		cout << endl << endl;
		cout << " Serial Dijkstra time: " << t2_step - t1_step << endl;
		if (test_for_iter_algoritm(count_vertex)) {
			if (is_correct_implementation(d_step, d_parallel, count_vertex)) {
				cout << "parallel implementayion is correct" << endl;
			}
			else {
				cout << "parallel implementayion isn`t correct" << endl;
			}
		}*/
	}
	MPI_Finalize();
}