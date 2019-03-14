#include "network_sim_parallel.h"
#include <sstream>

const int TASK_LEN = 6;

network_simiulation_parallel::network_simiulation_parallel() :
	corecount(0),
	rank(0),
	network_size(0),
	realizations(10),
	start_time(-1)
{}

void network_simiulation_parallel::MPI_main()
{
#ifdef _MPI
	start_time = MPI_Wtime();
#endif
	if (rank == 0)
		controlProcess();
	else if (rank > 0)
		computingProcess();
}

void network_simiulation_parallel::controlProcess()
{
#ifdef _MPI
	MPI_Status status, cur_status;
	double mpi_start_time = MPI_Wtime();

	cout << endl << "control_process() started" << endl;
	cout << "corecount " << corecount << endl;
	cout << "rank " << rank << endl;
	cout << "network_size " << network_size << endl;
	cout << "realizations " << realizations << endl;
	
	N_vec.push_back(network_size);
	c_vec.push_back(0.05);
	c_vec.push_back(0.1);
	c_vec.push_back(0.25);
	c_vec.push_back(0.5);
	k_vec.push_back(4);
	k_vec.push_back(6);
	k_vec.push_back(8);
	q_vec.push_back(1);
	for (int i = 0; i < 20; i++) // p = {0.0, 0.05, 0.1, ..., 0.95, 1.0}
		p_vec.push_back((double)i*0.05);
	for (int i = 0; i < realizations; i++)
		realization_vec.push_back(i);
	
	vector<vector<double>> search_space_params(6);
	for (auto x : q_vec)
		search_space_params[0].push_back((double)x);
	for (auto x : k_vec)
		search_space_params[1].push_back(x);
	for (auto x : c_vec)
		search_space_params[2].push_back(x);
	for (auto x : N_vec)
		search_space_params[3].push_back((double)x);
	for (auto x : p_vec)
		search_space_params[4].push_back(x);
	for (auto x : realization_vec)
		search_space_params[5].push_back((double)x);
	
	cout << "search_space_params : \n";
	for (auto x : search_space_params) {
		for (auto y : x)
			cout << y << " ";
		cout << "\n";
	}

	vector<point> search_space_points = generateSearchSpace(search_space_params);
	vector<double> tasks_times(search_space_points.size());
	for (auto x : tasks_times)
		x = -1.0;
	cout << "search space size : " << search_space_points.size() << endl;

	cout << "corecount " << corecount << endl;
	cout << "first " << corecount << " points from the search space :\n";
	for (int i = 0; i < corecount; i++) {
		for (auto y : search_space_points[i])
			cout << y << " ";
		cout << "\n";
	}
	
	double *task = new double[TASK_LEN];
	int task_index_for_sending = 0;

	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		for (unsigned j = 0; j < search_space_points[task_index_for_sending].size(); j++)
			task[j] = search_space_points[task_index_for_sending][j];
		MPI_Send(&task_index_for_sending, 1, MPI_INT, computing_process_index, 0, MPI_COMM_WORLD);
		MPI_Send(task, TASK_LEN, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
		tasks_times[task_index_for_sending] = MPI_Wtime();
		task_index_for_sending++;
	}
	cout << task_index_for_sending << " tasks were sent" << endl;

	int processed_task_count = 0;
	int received_task_index = -1;
	int stop_message = -1;
	string tasks_times_file_name = "tasks_times_out";
	
	while (processed_task_count < search_space_points.size()) {
		MPI_Recv(&received_task_index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		processed_task_count++;
		tasks_times[received_task_index] = MPI_Wtime() - tasks_times[received_task_index];

		ofstream ofile(tasks_times_file_name, ios_base::out);
		for (unsigned i=0; i<search_space_points.size(); i++) {
			for (unsigned j = 0; j < search_space_points[i].size(); j++) {
				ofile << search_space_points[i][j];
				if (j < search_space_points[i].size() - 1)
					ofile << "_";
			}
			if ( (tasks_times[i] < 1e6) && (tasks_times[i] > 0.0) )
				ofile << " " << tasks_times[i] << " seconds";
			else
				ofile << " -1";
			ofile << "\n";
		}
		ofile.close();
		
		cout << "processed_task_count " << processed_task_count;
		cout << " , time from start " << MPI_Wtime() - mpi_start_time << " s" << endl;
		
		if (task_index_for_sending >= search_space_points.size())
			MPI_Send(&stop_message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		else {
			for (unsigned j = 0; j < search_space_points[task_index_for_sending].size(); j++)
				task[j] = search_space_points[task_index_for_sending][j];
			MPI_Send(&task_index_for_sending, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			MPI_Send(task, TASK_LEN, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			tasks_times[task_index_for_sending] = MPI_Wtime();
			task_index_for_sending++;
		}
	}
	cout << "final time : " << MPI_Wtime() - mpi_start_time << endl;
	cout << "finilizing process " << rank << endl;
	MPI_Finalize();
#endif
}

vector<point> network_simiulation_parallel::generateSearchSpace(vector<vector<double>> search_space_params)
{
	vector<point> search_space_points;
	vector<int> index_arr;
	vector<double> cur_point;
	while (next_cartesian(search_space_params, index_arr, cur_point))
		search_space_points.push_back(cur_point);
	return search_space_points;
}

void network_simiulation_parallel::computingProcess()
{
	cout << "computingProcess()\n";
	double *task = new double[TASK_LEN];
	int task_index = -1;
#ifdef _MPI
	MPI_Status status;
	for (;;) {
		MPI_Recv(&task_index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if (task_index == -1) {
			cout << "finilizing process " << rank << endl;
			MPI_Finalize();
			break;
		}
		MPI_Recv(task, TASK_LEN, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		/*cout << "received task index " << task_index << endl;
		cout << "received task: \n";
		for (unsigned i = 0; i < TASK_LEN; i++)
			cout << task[i] << " ";
		cout << endl;*/

		network_simiulation_sequential n_s_seq;
		n_s_seq.verbosity = 0;
		n_s_seq.model = "qvoter_same";
		n_s_seq.network_type = "er";
		n_s_seq.seed = (unsigned int)(time(NULL)) * task_index;
		// "_q-" << q << "_k-" << k << "_c-" << c << "_N-" << N << "_p-" << p;
		n_s_seq.q = (int)task[0];
		n_s_seq.k = task[1];
		n_s_seq.c = task[2];
		n_s_seq.N = (int)task[3];
		n_s_seq.p = task[4];
		n_s_seq.realization = (int)task[5];

		n_s_seq.GetOutputName();
		n_s_seq.Init();
		n_s_seq.CreateGraphER();
		n_s_seq.CreateNodesState();
		n_s_seq.LaunchSimulation();
		n_s_seq.SaveMeasure();

		MPI_Send(&task_index, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
	}

	delete[] task;

	cout << "finilizing process " << rank << endl;
	MPI_Finalize();
#endif
}
