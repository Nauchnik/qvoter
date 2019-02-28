#ifndef network_sim_parallel_h
#define network_sim_parallel_h

#ifdef _MPI
#include <mpi.h>
#endif

#include "network_sim.h"

class network_simiulation_parallel: public network_simiulation_sequential
{
public:
	int corecount;
	int rank;
	double mpi_start_time;
	void MPI_main();
private:
	void controlProcess();
	void computingProcess();
};

#endif