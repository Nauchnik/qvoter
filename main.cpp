/* 
qvoter dynamics simulation

authors: 
Joanna Toruniewska, Warsaw University of Technology
Oleg Zaikin, ITMO University 
*/

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <cstdlib>
#include <math.h>
#include <chrono>

#include "network_sim.h"
#include "network_sim_parallel.h"

int main(int argc, char **argv) 
{
#ifdef _MPI
	int rank = 0;
	int corecount = 1;

	if (argc < 3) {
		cerr << "Usage: prog network-size realizations\n";
		exit(-1);
	}

	int network_size = atoi(argv[1]);
	int realizations = atoi(argv[2]);
	
	// parallel mode
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &corecount);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	network_simiulation_parallel n_s_par;
	n_s_par.verbosity = 0;
	n_s_par.network_size = network_size;
	n_s_par.realizations = realizations;
	n_s_par.corecount = corecount;
	n_s_par.rank = rank;
	n_s_par.MPI_main();
#else

	network_simiulation_sequential n_s_seq;
#ifdef _DEBUG
	n_s_seq.verbosity = 2;
	n_s_seq.model = "qvoter_same";
	n_s_seq.network_type = "er";
	// "_q-" << q << "_k-" << k << "_c-" << c << "_N-" << N << "_p-" << p;
	n_s_seq.q = 1;
	n_s_seq.k = 4;
	n_s_seq.c = 0.1;
	n_s_seq.N = 100000;
	n_s_seq.p = 0.25;
	n_s_seq.realization = 2;
#else
	if (argc < 12) {
		cerr << "Usage: program model network N c k q p t_0 t_max folder filename [realization] [seed]" << "\n";
		cerr << "--- few minutes execution ---\n"
			<< "model = qvoter_same\n"
			<< "network = er\n"
			<< "N = 400\n"
			<< "c = 0.5\n"
			<< "k = 8\n"
			<< "q = 1\n"
			<< "p = 0.55\n"
			<< "t0 = 0\n"
			<< "t_max = 100000\n"
			<< "folder = ./\n"
			<< "filename = 1\n";
		exit(-1);
	}
	n_s_seq.ReadParams(argc, argv);
#endif
	auto start = chrono::system_clock::now();

	n_s_seq.Init();
	n_s_seq.LaunchSimulation();
	n_s_seq.SaveMeasure();

	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "elapsed time: " << elapsed_seconds.count() << " s\n";
#endif

	return 0;
}

