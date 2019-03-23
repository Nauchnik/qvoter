#include "network_sim_parallel.h"
#include <sstream>

const int TASK_LEN = 6;

network_simiulation_parallel::network_simiulation_parallel() :
	corecount(0),
	rank(0),
	verbosity(1),
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

	vector<space_point> search_space_points = generateSearchSpace(search_space_params);
	cout << "search space size : " << search_space_points.size() << endl;

	string tasks_times_file_name = "tasks_times_out";
	ifstream tasks_file(tasks_times_file_name);
	
	vector<task> tasks_vec;
	if (tasks_file.is_open()) {
		cout << "opened tasks file " << tasks_times_file_name << endl;
		string str, word1, word2;
		task cur_task;
		while (getline(tasks_file, str)) {
			if (str.size() < 2)
				continue;
			stringstream sstream;
			sstream << str;
			sstream >> word1 >> word2;
			cur_task.solving_time = -1.0;
			if ( (word2 == "-1") || (word2 == "0") )
				cur_task.status = -1; // not launched
			else if (word2 == "-2")
				cur_task.status = -2; // interrupted
			else {
				cur_task.status = 1; // solved
				istringstream(word2) >> cur_task.solving_time;
			}
			cur_task.out_file_name = "";
			tasks_vec.push_back(cur_task);
		}
		cout << "tasks_vec.size() " << tasks_vec.size() << endl;
		cout << "first 20 unsolved tasks statuses and times: \n";
		int uns = 0;
		for (auto task : tasks_vec) {
			if (task.status != 1) {
				cout << task.status << " " << task.solving_time << endl;
				uns++;
				if (uns == 20)
					break;
			}
		}
	}
	else {
		tasks_vec.resize(search_space_points.size());
		for (auto x : tasks_vec)
			x.status = -1;
	}
	tasks_file.close();

	if (search_space_points.size() != tasks_vec.size()) {
		cerr << "search_space_points.size() != tasks_vec.size()" << endl;
		cerr << search_space_points.size() << " != " << tasks_vec.size() << endl;
		MPI_Abort(MPI_COMM_WORLD, 0);
	}
	for (unsigned i = 0; i < search_space_points.size(); i++)
		tasks_vec[i].task_point = search_space_points[i];
	search_space_points.clear();
	
	cout << "corecount " << corecount << endl;
	
	double *task = new double[TASK_LEN];
	int task_index_for_sending = 0;
	int skipped_tasks = 0;

	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		while (tasks_vec[task_index_for_sending].status == 1) {
			task_index_for_sending++;
			skipped_tasks++;
		}
		for (unsigned j = 0; j < tasks_vec[task_index_for_sending].task_point.size(); j++)
			task[j] = tasks_vec[task_index_for_sending].task_point[j];
		MPI_Send(&task_index_for_sending, 1, MPI_INT, computing_process_index, 0, MPI_COMM_WORLD);
		MPI_Send(task, TASK_LEN, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
		tasks_vec[task_index_for_sending].solving_time = MPI_Wtime();
		task_index_for_sending++;
	}
	cout << skipped_tasks << " tasks were skipped" << endl;
	cout << task_index_for_sending << " tasks were sent" << endl;
	
	int processed_task_count = 0;
	int received_task_index = -1;
	int stop_message = -1;
	int task_status = -1;
	int previous_launches_count = 0;
	while (processed_task_count + skipped_tasks < tasks_vec.size()) {
		MPI_Recv(&received_task_index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&task_status, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&previous_launches_count, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		processed_task_count++;
		tasks_vec[received_task_index].solving_time = previous_launches_count*MAX_SOLVING_TIME_SEC + 
			(MPI_Wtime() - tasks_vec[received_task_index].solving_time);
		tasks_vec[received_task_index].status = task_status;
		cout << "task_index " << received_task_index << endl;
		cout << "solving_time " << tasks_vec[received_task_index].solving_time << endl;
		cout << "task_status " << task_status << endl;
		
		ofstream ofile(tasks_times_file_name, ios_base::out);
		for (unsigned i=0; i< tasks_vec.size(); i++) {
			for (unsigned j = 0; j < tasks_vec[i].task_point.size(); j++) {
				ofile << tasks_vec[i].task_point[j];
				if (j < tasks_vec[i].task_point.size() - 1)
					ofile << "_";
			}
			//if ( (tasks_vec[i].solving_time < 1e6) && (tasks_vec[i].solving_time > 0.0) )
			if ((tasks_vec[i].status == 1) && (tasks_vec[i].solving_time < 1e6) && (tasks_vec[i].solving_time > 0.0))
				ofile << " " << tasks_vec[i].solving_time << " seconds";
			else
				ofile << " " << tasks_vec[i].status;
			//else
			//	ofile << " -1";
			ofile << "\n";
		}
		ofile.close();
		
		cout << "processed_task_count=" << processed_task_count;
		cout << " , skipped_tasks=" << skipped_tasks;
		cout << " , time from start " << MPI_Wtime() - mpi_start_time << " s" << endl;
		
		while (tasks_vec[task_index_for_sending].status == 1) {
			task_index_for_sending++;
			skipped_tasks++;
		}
		
		if (task_index_for_sending >= tasks_vec.size())
			MPI_Send(&stop_message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		else {
			for (unsigned j = 0; j < tasks_vec[task_index_for_sending].task_point.size(); j++)
				task[j] = tasks_vec[task_index_for_sending].task_point[j];
			MPI_Send(&task_index_for_sending, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			MPI_Send(task, TASK_LEN, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			tasks_vec[task_index_for_sending].solving_time = MPI_Wtime();
			task_index_for_sending++;
		}
	}
	cout << "final time : " << MPI_Wtime() - mpi_start_time << endl;
	cout << "finilizing process " << rank << endl;
	MPI_Finalize();
#endif
}

vector<space_point> network_simiulation_parallel::generateSearchSpace(vector<vector<double>> search_space_params)
{
	vector<space_point> search_space_points;
	vector<int> index_arr;
	vector<double> cur_point;
	while (next_cartesian(search_space_params, index_arr, cur_point))
		search_space_points.push_back(cur_point);
	return search_space_points;
}

void network_simiulation_parallel::computingProcess()
{
	//cout << "computingProcess()\n";
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
		if (verbosity > 1) {
			cout << "received task index " << task_index << endl;
			cout << "received task: \n";
			for (unsigned i = 0; i < TASK_LEN; i++)
				cout << task[i] << " ";
			cout << endl;
		}
		
		network_simiulation_sequential n_s_seq;
		n_s_seq.verbosity = 0;
		n_s_seq.model = "qvoter_same";
		n_s_seq.network_type = "er";
		// "_q-" << q << "_k-" << k << "_c-" << c << "_N-" << N << "_p-" << p;
		n_s_seq.q = (int)task[0];
		n_s_seq.k = task[1];
		n_s_seq.c = task[2];
		n_s_seq.N = (int)task[3];
		n_s_seq.p = task[4];
		n_s_seq.realization = (int)task[5];

		string out_f_name_wout_seed = n_s_seq.GetOutputNameWoutSeed();
		string step_file_name = "";
		unsigned seed = 0;
		n_s_seq.FindStateFileName(out_f_name_wout_seed, step_file_name, seed);
		if (verbosity > 1) {
			cout << "step_file_name " << step_file_name << endl;
			cout << "seed " << seed << endl;
		}
		
		if ((step_file_name == "") && (seed == 0))
			n_s_seq.seed = (unsigned int)(time(NULL)) * task_index;
		else
			n_s_seq.seed = seed;
		
		n_s_seq.GetOutputName();
		n_s_seq.start_time = MPI_Wtime();
		n_s_seq.Init();
		if (step_file_name != "") { // read state from file
			cout << "Trying to read file " << step_file_name << endl;
			n_s_seq.ReadSimulationState(step_file_name);
		}
		else { // start from scratch
			n_s_seq.CreateGraphER();
			n_s_seq.CreateNodesState();
		}
		n_s_seq.LaunchSimulation();
		if (n_s_seq.status != -1) // if solved or interrupted
			n_s_seq.SaveMeasure();
		
		MPI_Send(&task_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&n_s_seq.status, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&n_s_seq.previous_launches_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	delete[] task;

	cout << "finilizing process " << rank << endl;
	MPI_Finalize();
#endif
}
