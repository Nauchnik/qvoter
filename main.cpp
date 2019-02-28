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

int main(int argc, char **argv) 
{
	if (argc < 11) {
		cerr << "Usage: program model network N c k q p t_max folder filename [seed]" << "\n";
		cerr << "--- few minutes execution ---\n"
			<< "model = qvoter_same\n"
			<< "network = er\n"
			<< "N = 400\n"
			<< "c = 0.5\n"
			<< "k = 8\n"
			<< "q = 1\n"
			<< "p = 0.55\n"
			<< "t_max = 100000\n"
			<< "folder = ./\n"
			<< "filename = 1\n";
		exit(-1);
	}

	auto start = std::chrono::system_clock::now();

	network_simiulation_sequential n_s_s;
	n_s_s.ReadParams(argc, argv);
	n_s_s.Init();
	n_s_s.CreateGraphER();
	n_s_s.CreateNodesState();
	n_s_s.LaunchSimulation();
	n_s_s.SaveMeasure();

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";

	return 0;
}

