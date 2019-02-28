#include "network_sim_parallel.h"
#include <sstream>

const int TASK_LEN = 10;

void network_simiulation_parallel::MPI_main()
{
	if (rank == 0)
		controlProcess();
	else if (rank > 0)
		computingProcess();
}

void network_simiulation_parallel::controlProcess()
{
#ifdef _MPI
	MPI_Status status, cur_status;
	stringstream sstream_out;
	mpi_start_time = MPI_Wtime();

	cout << endl << "control_process() started" << endl;

	double *task = new double[TASK_LEN];

	int send_task_count = 0;

	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		sendTaskIls(task, send_task_count, computing_process_index, depths_vec[send_task_count]);
		send_task_count++;
	}
	sstream_out << "send_task_count " << send_task_count << endl;
	sstream_out << "tasks count " << depths_vec.size() << endl;

	writeOutputData(sstream_out);

	unsigned processed_task_count = 0;
	int received_task_index;
	double received_residual;
	double task_processing_time;
	int stop_message = -1;

	// get results and send new tasks on idle computing processes
	while (processed_task_count < depths_vec.size()) {
		MPI_Recv(&received_task_index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		cur_status = status;
		MPI_Recv(&task_processing_time, 1, MPI_DOUBLE, cur_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(result, ILS_RESULT_LEN, MPI_DOUBLE, cur_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		processed_task_count++;

		sstream_out << "processed_task_count " << processed_task_count;
		sstream_out << " , time from start " << MPI_Wtime() - mpi_start_time << " s" << endl;

		received_residual = result[0];

		if (received_residual < record_point.residual) {
			record_point.residual = received_residual;
			record_point.cb = result[1];
			record_point.rhob = result[2];
			record_point.R = result[3];
			record_point.tau = result[4];
			unsigned cur_depths_size = depths_vec[received_task_index].size();
			record_point.depths.resize(cur_depths_size);
			unsigned cur_cws_size = cur_depths_size - 1;
			record_point.cws.resize(cur_cws_size);
			for (unsigned j = 0; j < cur_cws_size; j++)
				record_point.cws[j] = result[5 + j];
			for (unsigned j = 0; j < cur_depths_size; j++)
				record_point.depths[j] = result[5 + cur_cws_size + j];

			if (record_point.depths != depths_vec[received_task_index]) {
				cerr << "record_point.depths != depths_vec[received_task_index]" << endl;
				cerr << "received_task_index " << received_task_index << endl;
				MPI_Abort(MPI_COMM_WORLD, 0);
				return;
			}

			//cout << "Control process, new residual minimum : "  << received_residual << endl;
			sstream_out << endl << "Control process, new residual minimum:" << endl;
			sstream_out << "task processing time " << task_processing_time << " s" << endl;
			sstream_out << "task index " << received_task_index << endl;
			sstream_out << strPointData(record_point);
			sstream_out << "time from start " << MPI_Wtime() - mpi_start_time << " s" << endl;
		}
		// if free tasks for sending
		if (send_task_count < depths_vec.size()) {
			sendTaskIls(task, send_task_count, status.MPI_SOURCE, depths_vec[send_task_count]);
			send_task_count++;
		}
		else
			MPI_Send(&stop_message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);

		writeOutputData(sstream_out);
	}

	sstream_out << endl << "SEARCH ENDED!" << endl;
	sstream_out << strPointData(record_point);
	sstream_out << "final time " << MPI_Wtime() - mpi_start_time << " s" << endl;

	writeOutputData(sstream_out);

	delete[] task;
	delete[] result;

	cout << "finilizing process " << rank << endl;
	MPI_Finalize();
#endif
}

void network_simiulation_parallel::computingProcess()
{
#ifdef _MPI

#endif
}
