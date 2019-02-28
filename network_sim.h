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
	vector<vector<double>> search_space_values;
	vector<vector<double>> search_space_points;

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
	void GenerateSearchSpace();
};

template< typename T >
bool next_cartesian(vector<T> &vii, vector<int> &index_arr, T &cur_vi)
{
	if (index_arr.size() == 0) { // init
		index_arr.resize(vii.size());
		//for( auto &x : index_arr )
		//	x = 0;
		for (vector<int> ::iterator it = index_arr.begin(); it != index_arr.end(); ++it)
			*it = 0;
	}
	if (index_arr[0] == -1)
		return false;
	// get current value
	cur_vi.resize(vii.size());
	for (unsigned i = 0; i < index_arr.size(); ++i)
		cur_vi[i] = vii[i][index_arr[i]];
	// check if last iteration
	bool IsLastValue = true;
	for (unsigned i = 0; i < index_arr.size(); ++i) {
		if (index_arr[i] != vii[i].size() - 1) {
			IsLastValue = false;
			break;
		}
	}
	if (IsLastValue)
		index_arr[0] = -1; // condition of stopping
	else {
		// find last changable row to increase its value
		int last_changable = index_arr.size() - 1;
		while (last_changable != -1) {
			if (index_arr[last_changable] < (int)(vii[last_changable].size() - 1))
				break;
			--last_changable;
		}
		index_arr[last_changable]++;
		for (unsigned i = last_changable + 1; i < index_arr.size(); ++i)
			index_arr[i] = 0;
	}

	return true;
}

#endif