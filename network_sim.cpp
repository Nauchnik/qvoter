#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <vector>
#include <list>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <string>

using namespace std;

//creates random Erdos Renyi graphwith N nodes and mean degree k
void CreateGraphER(vector < vector <int> > &network, int N, double k){
		double p=k/(N-1);
		double r;
		for (int i=0; i<N;i++){
			for(int j=i; j<N;j++){
				r =((double)rand()) / RAND_MAX;
				if(r<p){
					network[i].push_back(j);
					network[j].push_back(i);
				}
			}
		}
}
//adds states to nodes, c*N nodes with plus state, (1-c)*N with -1 state
int CreateNodesState(int N, double c, int*nodes_states){
	for (int i=0; i<round(c*N); i++){
		nodes_states[i]=1;
	}
	for (int i=round(c*N); i<N; i++){
		nodes_states[i]=-1;
	}
	return 0;
}

struct single_measure{
	int E_interface;
	int E_plus;
	int E_minus;
	int plus_nodes;
	int minus_nodes;
	int time;
};

//measure number of links in each state and number of nodes in each state
single_measure Measure(int* nodes_states, int N, vector < vector <int> > &network, int time){ //przemyslec co z tym wskaznikiem, int* E_int){
	double plus_nodes=0;
	double minus_nodes=0;
	int E_int=0;
	int E_plus=0;
	int E_minus=0;
	for (int i=0; i<N; i++){
		int neighbor_number=network[i].size();
		if(nodes_states[i] ==1){
			plus_nodes++;
			for(int j=0;j<neighbor_number;j++){
				if(nodes_states[network[i][j]]==1){
					E_plus++;
				}else{
					E_int++;
				}
			}
		}else{
				minus_nodes++;
				for(int j=0;j<neighbor_number;j++){
				if(nodes_states[network[i][j]]==-1){
					E_minus++;
				}else{
					E_int++;
				}
			}
		}
	}
	single_measure cokolwiek;
	cokolwiek.plus_nodes=plus_nodes;
	cokolwiek.minus_nodes=minus_nodes;
	cokolwiek.E_interface=E_int/2;
	cokolwiek.E_plus=E_plus/2;
	cokolwiek.E_minus=E_minus/2;
	cokolwiek.time=time;

	return cokolwiek;
}

//check if two nodes are connected
bool isConnected(vector < vector <int> > &network,int node1, int node2){
	bool connected=false;
	if (std::find(network[node1].begin(), network[node1].end(),node2)!=network[node1].end()){
		connected=true;
	}
	return connected;
}
//rewire basic_node link from old_neighbor to randomly chosen node in the network
void RewireToRandom(vector < vector <int> > &network, int basic_node,int old_neighbor, int N ){
	//TODO lack of ensuracne that there is no more not connected nodes in the network
	int new_neighbor=rand() %N;
	while(isConnected(network, basic_node, new_neighbor)){
		new_neighbor=rand() %N;
	}
	replace (network[basic_node].begin(), network[basic_node].end(),old_neighbor , new_neighbor);
	std::vector<int>::iterator position = std::find(network[old_neighbor].begin(), network[old_neighbor].end(), basic_node);
	network[old_neighbor].erase(position);
	network[new_neighbor].push_back(basic_node);

}

//rewire basic_node link from old_neighbor to randomly chosen node in the network with the same state as basic_node
void RewireToSame(vector < vector <int> > &network, int basic_node,int old_neighbor, int* nodes_states, int inTheSameState){
//very very slow probably
	int neighbor_in_the_same_state=0;
	int state_basic_nodes=nodes_states[basic_node];
	for(int i=0;i<network[basic_node].size();i++){
		if(nodes_states[network[basic_node][i]]==state_basic_nodes){
			neighbor_in_the_same_state++;
		}
	}
	int number_of_possible_nodes=inTheSameState-neighbor_in_the_same_state;
	if(number_of_possible_nodes>0){
		int rand_number=rand() %(number_of_possible_nodes);
		int j=0, new_neighbor=-1;
		while(j<=rand_number){
			new_neighbor++;
			if(nodes_states[new_neighbor]==state_basic_nodes && !isConnected(network, basic_node, new_neighbor)){
				j++;
			}
		}

		replace (network[basic_node].begin(), network[basic_node].end(),old_neighbor , new_neighbor);
		std::vector<int>::iterator position = std::find(network[old_neighbor].begin(), network[old_neighbor].end(), basic_node);
		network[old_neighbor].erase(position);
		network[new_neighbor].push_back(basic_node);
	}
}



bool check_qpanel(vector<int> &qpanel, int* nodes_states, int q){
	bool identical=true;
	int state=nodes_states[qpanel[0]];
	for(int i=1; i<q;i++){
		if(nodes_states[qpanel[i]]!=state){
			identical=false;
			break;
		}
	}
	return identical;
}
//whole random qvoter model dynamics
void DynamicsQVoterRandom(list<single_measure> &measure_list, int t_max,
		vector<vector<int> > &network, int* nodes_states, int N, int q,
		double p) {

	int t = 0;
	vector<int> qpanel(q);
	int basic_node, r1, old_neighbor;

	int size;
	measure_list.push_back(Measure(nodes_states, N, network, t)); //, &E_int));
	while ((t < t_max) && (measure_list.back().E_interface > 0)) { //E_int>0){
		for (int microstep = 0; microstep < N; microstep++) {
			basic_node = rand() % N;
			//TODO zliczanie?
			size = network[basic_node].size();
			if (size > 0) {
				r1 = rand() % size;
				qpanel[0] = network[basic_node][r1];
				if (nodes_states[basic_node] != nodes_states[qpanel[0]]) {
					for (int i = 1; i < q; i++) {
						r1 = rand() % size;
						qpanel[i] = network[basic_node][r1];
					};
					if (check_qpanel(qpanel, nodes_states, q)) {
						if ((((double) rand()) / RAND_MAX) < p) {
							r1 = rand() % q;
							old_neighbor = qpanel[r1];

							RewireToRandom(network, basic_node, old_neighbor,
									N);
						} else {
							nodes_states[basic_node] = nodes_states[qpanel[0]];
						}
					}
				}
			}
		}
		t++;
		measure_list.push_back(Measure(nodes_states, N, network, t)); //, &E_int));
		if ((t % 100) == 0)
			cout << t << "\n";
	}

}

//whole same qvoter model dynamics
void DynamicsQVoterSame(list<single_measure> &measure_list, int t_max,
		vector<vector<int> > &network, int* nodes_states, int N, int q,
		double p) {
	int t = 0;
	vector<int> qpanel(q);
	int basic_node, r1, old_neighbor, inTheSameState;
	int number_of_plus_nodes = 0, number_of_minus_nodes = 0;
	for (int i = 0; i < N; i++) {
		if (nodes_states[i] == 1) {
			number_of_plus_nodes++;
		} else {
			number_of_minus_nodes++;
		}
	}
	int size;
	measure_list.push_back(Measure(nodes_states, N, network, t)); //, &E_int));
	while ((t < t_max) && (measure_list.back().E_interface > 0)) { //E_int>0){
		for (int microstep = 0; microstep < N; microstep++) {
			basic_node = rand() % N;
			//TODO zliczanie?
			size = network[basic_node].size();
			if (size > 0) {
				r1 = rand() % size;
				qpanel[0] = network[basic_node][r1];
				if (nodes_states[basic_node] != nodes_states[qpanel[0]]) {
					for (int i = 1; i < q; i++) {
						r1 = rand() % size;
						qpanel[i] = network[basic_node][r1];
					};
					if (check_qpanel(qpanel, nodes_states, q)) {
						if ((((double) rand()) / RAND_MAX) < p) {
							r1 = rand() % q;
							old_neighbor = qpanel[r1];

							if (nodes_states[basic_node] == 1) {
								inTheSameState = number_of_plus_nodes;
							} else {
								inTheSameState = number_of_minus_nodes;
							}
							RewireToSame(network, basic_node, old_neighbor,
									nodes_states, inTheSameState);
						} else {
							if (nodes_states[basic_node] == 1) {
								number_of_plus_nodes--;
								number_of_minus_nodes++;
							} else {
								number_of_plus_nodes++;
								number_of_minus_nodes--;
							}
							nodes_states[basic_node] = nodes_states[qpanel[0]];
						}
					}
				}
			}
		}
		t++;
		measure_list.push_back(Measure(nodes_states, N, network, t)); //, &E_int));
		if ((t % 100) == 0)
			cout << t << "\n";
	}
}
//write measure list to the file
void SaveMeasure( const std::list<single_measure>& measure_list, string  outfile_name){
	{
	    std::ofstream file(outfile_name);
	    file << "#time\tE_interface\tE_plus\tE_minus\tplus_nodes\tminus_nodes\n";
	    for( const auto& v : measure_list ) file <<v.time<<'\t'<< v.E_interface << '\t'<<v.E_plus<<'\t'<<v.E_minus<<'\t'<<v.plus_nodes<<'\t'<<v.minus_nodes<<'\n' ;
	}
}

void create_params() {
	string catalogue = "";
	string name = "params_";

	for (int i = 0; i < 10; i++) {
		std::ofstream file(catalogue + name + to_string(i));
		file << "model=qvoter_same" << '\n';
		file << "network=er" << '\n';
		file << "N=1000" << '\n';
		file << "c=0.5" << '\n';
		file << "k=4" << '\n';
		file << "q=1" << '\n';
		file << "t_max=20000000" << '\n';
		file << "p=0.0" << '\n';
		file << "folder=/home/joanna/workspace-cdt/network_sim/" << '\n';
		file << "filename=test_" + to_string(i) << '\n';
	}
}


int main(int argc, char **argv)
{


	if (argc > 1) {
		string params = argv[1];

		std::srand(std::time(0));

		int N;
		double c;
		double k;
		double p;
		int q;
		int t_max;
		string model;
		std::string filename;
		std::string folder;
		string outname;
		list<single_measure> measure_list = { };

// params file should have:
// value = param

		std::ifstream cFile(params);
		cout<<params<<"\t";
		if (cFile.is_open()) {
			std::string line;
			while (getline(cFile, line)) {
				//line.erase(remove_if(line.begin(), line.end(), isspace),
				//                   line.end());
				if (line[0] == '#' || line.empty())
					continue;
				auto delimiterPos = line.find("=");
				auto name = line.substr(0, delimiterPos);
				auto value = line.substr(delimiterPos + 1);
				std::cout << "#" << name << " " << value << '\n';
				if (name == "model") {
					model = (value);
				}
				if (name == "N") {
					N = std::stoi(value);
				}
				if (name == "c") {
					c = std::stod(value);
				}
				if (name == "k") {
					k = std::stod(value);
				}
				if (name == "p") {
					p = std::stod(value);
				}
				if (name == "q") {
					q = std::stoi(value);
				}
				if (name == "t_max") {
					t_max = std::stoi(value);
				}
				if (name == "filename") {
					filename = value;
				}
				if (name == "folder") {
					folder = value;
				}
			}

			outname = folder + model + "_q-" + to_string(q) + "_k-"
					+ to_string(k) + "_c-" + to_string(c) + "_N-" + to_string(N)
					+ "_p-" + to_string(p) + "_" + filename;

			int nodes_states[N];
			vector<vector<int> > network(N, vector<int>());

			CreateGraphER(network, N, k);
			CreateNodesState(N, c, nodes_states);

			if (model == "qvoter_same") {
				DynamicsQVoterSame(measure_list, t_max, network, nodes_states,
						N, q, p);
			} else if (model == "qvoter_random") {
				DynamicsQVoterRandom(measure_list, t_max, network, nodes_states,
						N, q, p);
			}
			SaveMeasure(measure_list, outname);
		} else {
			std::cerr << "Couldn't open config file for reading.\n";
		}
	} else {
		cout << "Enter name of the file with params" << "\n";
	}
	return 0;
}


