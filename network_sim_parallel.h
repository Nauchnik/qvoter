#ifndef network_sim_parallel_h
#define network_sim_parallel_h

#ifdef _MPI
#include <mpi.h>
#endif

#include "network_sim.h"

#define point vector<double>

class network_simiulation_parallel: public network_simiulation_sequential
{
public:
	network_simiulation_parallel();
	int corecount;
	int rank;
	int network_size;
	vector<int> q_vec;
	vector<double> k_vec;
	vector<double> c_vec;
	vector<int> N_vec;
	vector<double> p_vec;
	vector<int> realization_vec;

	void MPI_main();

private:
	void controlProcess();
	vector<point> generateSearchSpace(vector<vector<double>> search_space_values);
	void computingProcess();
};

#endif