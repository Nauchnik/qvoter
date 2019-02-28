#ifndef network_sim_h
#define network_sim_h

#include <vector>
#include <string>
#include <list>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

using namespace std;

struct single_measure {
	int E_interface;
	int E_plus;
	int E_minus;
	int plus_nodes;
	int minus_nodes;
	int time;
};

class network_simiulation_sequential
{
public:
	string params_file_name;
	unsigned int seed;
	int N;
	double c;
	double k;
	double p;
	int q;
	int t_max;
	string model;
	string filename;
	string folder;
	string outname;
	list<single_measure> measure_list;
	vector<int> nodes_states;
	vector<vector<int> > network;

	void ReadParams(const int argc, char **argv);
	void Init();
	void LaunchSimulation();
	void CreateGraphER();
	int CreateNodesState();
	void SaveMeasure();

protected:
	single_measure Measure(vector<int> &nodes_states, int N, vector < vector <int> > &network, int time);
	void RewireToRandom(vector < vector <int> > &network, int basic_node, int old_neighbor, int N);
	void RewireToSame(vector < vector <int> > &network, int basic_node, int old_neighbor,
		vector<int> &nodes_states, int inTheSameState);
	bool CheckQpanel(vector<int> &qpanel, vector<int> &nodes_states, int q);
	void GenerateParams();
	void CreateParams();
	bool isConnected(vector < vector <int> > &network, int node1, int node2);
	void DynamicsQVoterRandom(list<single_measure> &measure_list, int t_max,
		vector<vector<int> > &network, vector<int> &nodes_states, int N, int q,
		double p);
	void DynamicsQVoterSame(list<single_measure> &measure_list, int t_max,
		vector<vector<int> > &network, vector<int> &nodes_states, int N, int q,
		double p);
	void DynamicsQVoterSameAK(list<single_measure> &measure_list, int t_max,
		vector<vector<int> > &network, vector<int> nodes_states, int N, int q,
		double p);
};

#endif