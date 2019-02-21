/* 
qvoter dynamics simulation

authors: 
Joanna Toruniewska, Warsaw University of Technology
Oleg Zaikin, ITMO University 
*/

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <list>
#include <string>

using namespace std;

struct single_measure {
	int E_interface;
	int E_plus;
	int E_minus;
	int plus_nodes;
	int minus_nodes;
	int time;
};

void CreateGraphER(vector < vector <int> > &network, int N, double k);
int CreateNodesState(int N, double c, vector<int> &nodes_states);
single_measure Measure(vector<int> &nodes_states, int N, vector < vector <int> > &network, int time);
bool isConnected(vector < vector <int> > &network, int node1, int node2);
void RewireToRandom(vector < vector <int> > &network, int basic_node, int old_neighbor, int N);
void RewireToSame(vector < vector <int> > &network, int basic_node, int old_neighbor, vector<int> &nodes_states, int inTheSameState);
bool CheckQpanel(vector<int> &qpanel, vector<int> &nodes_states, int q);
void DynamicsQVoterRandom(list<single_measure> &measure_list, int t_max,
	vector<vector<int> > &network, vector<int> &nodes_states, int N, int q,
	double p);
void DynamicsQVoterSame(list<single_measure> &measure_list, int t_max,
	vector<vector<int> > &network, vector<int> &nodes_states, int N, int q,
	double p);
void SaveMeasure(const list<single_measure>& measure_list, string  outfile_name);
void GenerateParams();
void CreateParams();

int main(int argc, char **argv) {
	string params_file_name;
	unsigned int seed;
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

#ifdef _DEBUG
	params_file_name = "..//params-1-10-800";
	seed = 3; // seed
#else
	if (argc < 11) {
		cerr << "Usage: prog params [seed]" << "\n";
		exit(-1);
	}
	//parameters in the order: model network N c k q p t_max filename folder
	model = argv[1];
	N = atoi(argv[3]);
	c = atof(argv[4]);
	k = atof(argv[5]);
	q = atoi(argv[6]);
	p = atof(argv[7]);
	t_max = atoi(argv[8]);
	folder = argv[9];
	filename = argv[10];

	seed = (unsigned int) (time(NULL)); //does not work under my compilier wothout brakets
	if (argc > 11) {
		istringstream(argv[11]) >> seed;
		cout << "seed " << seed << endl;
	} else
		cout << "seed was chosen from the current time\n";
	srand(seed);
#endif
	cout << "seed " << seed << endl;
	cout << "params: model-" << model << "  N-" << N << "  c-" << c << "  k-"
			<< k << "  q-" << q << "  p-" << p << "  t_max-" << t_max
			<< "  filename-" << filename << "  folder-" << folder << "\n";

	list<single_measure> measure_list = { };
//	ifstream cFile(params_file_name);
//	if (!cFile.is_open()){
//		cerr << "Couldn't open config file " << params_file_name << " for reading.\n";
//		exit(-1);
//	}
//
//	string line;
//	while (getline(cFile, line)) {
//		//line.erase(remove_if(line.begin(), line.end(), isspace),
//		//                   line.end());
//		if (line[0] == '#' || line.empty())
//			continue;
//		auto delimiterPos = line.find("=");
//		auto name = line.substr(0, delimiterPos);
//		auto value = line.substr(delimiterPos + 1);
//		cout << "#" << name << " " << value << '\n';
//		if (name == "model")
//			model = (value);
//		if (name == "N")
//			N = stoi(value);
//		if (name == "c")
//			c = stod(value);
//		if (name == "k")
//			k = stod(value);
//		if (name == "p")
//			p = stod(value);
//		if (name == "q")
//			q = stoi(value);
//		if (name == "t_max")
//			t_max = stoi(value);
//		if (name == "filename")
//			filename = value;
//		if (name == "folder")
//			folder = value;
//	}
//	cFile.close();

	outname = folder + model + "_q-" + to_string(q) + "_k-" + to_string(k)
			+ "_c-" + to_string(c) + "_N-" + to_string(N) + "_p-" + to_string(p)
			+ "_" + filename;
	cout << "output file name is " << outname << endl;

	vector<int> nodes_states(N);
	vector<vector<int> > network(N, vector<int>());

	CreateGraphER(network, N, k);
	CreateNodesState(N, c, nodes_states);

	if (model == "qvoter_same") {
		DynamicsQVoterSame(measure_list, t_max, network, nodes_states, N, q, p);
	} else if (model == "qvoter_random") {
		DynamicsQVoterRandom(measure_list, t_max, network, nodes_states, N, q,
				p);
	} else {
		cerr << "Wrong model name. Possible models: qvoter_same, qvoter_random"
				<< "\n";
		exit(-1);
	}

	SaveMeasure(measure_list, outname);

	return 0;
}

// create a random Erdos Renyi graph with N nodes and mean degree k
void CreateGraphER(vector < vector <int> > &network, int N, double k)
{
	double p=k/((double)N-1);
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
int CreateNodesState(int N, double c, vector<int> &nodes_states)
{
	int bound = (int)round(c*(double)N);
	for (int i=0; i<bound; i++)
		nodes_states[i]=1;
	for (int i=bound; i<N; i++)
		nodes_states[i]=-1;
	return 0;
}

//measure number of links in each state and number of nodes in each state
single_measure Measure(vector<int> &nodes_states, int N, vector < vector <int> > &network, int time)
{ //przemyslec co z tym wskaznikiem, int* E_int){
	int plus_nodes=0;
	int minus_nodes=0;
	int E_int=0;
	int E_plus=0;
	int E_minus=0;
	for (int i=0; i<N; i++){
		size_t neighbor_number=network[i].size();
		if(nodes_states[i] == 1){
			plus_nodes++;
			for(int j=0;j<neighbor_number;j++){
				if(nodes_states[network[i][j]]==1)
					E_plus++;
				else
					E_int++;
			}
		}else{
				minus_nodes++;
				for(int j=0;j<neighbor_number;j++){
				if(nodes_states[network[i][j]]==-1)
					E_minus++;
				else
					E_int++;
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
bool isConnected(vector < vector <int> > &network,int node1, int node2)
{
	bool connected=false;
	if (find(network[node1].begin(), network[node1].end(),node2)!=network[node1].end())
		connected=true;
	return connected;
}

//rewire basic_node link from old_neighbor to randomly chosen node in the network
void RewireToRandom(vector < vector <int> > &network, int basic_node,int old_neighbor, int N )
{
	//TODO lack of ensuracne that there is no more not connected nodes in the network
	int new_neighbor=rand() % N;
	while(isConnected(network, basic_node, new_neighbor))
		new_neighbor=rand() % N;
	replace (network[basic_node].begin(), network[basic_node].end(),old_neighbor , new_neighbor);
	vector<int>::iterator position = find(network[old_neighbor].begin(), network[old_neighbor].end(), basic_node);
	network[old_neighbor].erase(position);
	network[new_neighbor].push_back(basic_node);
}

// rewire basic_node link from old_neighbor to randomly chosen node in the network with the same state as basic_node
void RewireToSame(vector < vector <int> > &network, int basic_node,int old_neighbor, vector<int> &nodes_states, int inTheSameState)
{
//very very slow probably
	int neighbor_in_the_same_state=0;
	int state_basic_nodes=nodes_states[basic_node];
	for(unsigned i=0;i<network[basic_node].size();i++){
		if(nodes_states[network[basic_node][i]]==state_basic_nodes)
			neighbor_in_the_same_state++;
	}
	int number_of_possible_nodes=inTheSameState-neighbor_in_the_same_state;
	if(number_of_possible_nodes>0){
		int rand_number=rand() % (number_of_possible_nodes);
		int j=0, new_neighbor=-1;
		while(j<=rand_number){
			new_neighbor++;
			if(nodes_states[new_neighbor]==state_basic_nodes && !isConnected(network, basic_node, new_neighbor))
				j++;
		}

		replace (network[basic_node].begin(), network[basic_node].end(), old_neighbor, new_neighbor);
		vector<int>::iterator position = find(network[old_neighbor].begin(), network[old_neighbor].end(), basic_node);
		network[old_neighbor].erase(position);
		network[new_neighbor].push_back(basic_node);
	}
}

bool CheckQpanel(vector<int> &qpanel, vector<int> &nodes_states, int q)
{
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
		vector<vector<int> > &network, vector<int> &nodes_states, int N, int q,
		double p) 
{
	int t = 0;
	vector<int> qpanel(q);
	int basic_node, r1, old_neighbor;

	size_t size;
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
					}
					if (CheckQpanel(qpanel, nodes_states, q)) {
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
		vector<vector<int> > &network, vector<int> &nodes_states, int N, int q,
		double p) 
{
	int t = 0;
	vector<int> qpanel(q);
	int basic_node, r1, old_neighbor, inTheSameState;
	int number_of_plus_nodes = 0, number_of_minus_nodes = 0;
	for (int i = 0; i < N; i++) {
		if (nodes_states[i] == 1)
			number_of_plus_nodes++;
		else
			number_of_minus_nodes++;
	}
	size_t size;
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
					}
					if (CheckQpanel(qpanel, nodes_states, q)) {
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

// write the measure list to the output file
void SaveMeasure( const list<single_measure>& measure_list, string outfile_name )
{
	ofstream ofile(outfile_name);
	ofile << "#time\tE_interface\tE_plus\tE_minus\tplus_nodes\tminus_nodes\n";
	for( const auto& v : measure_list ) 
		ofile <<v.time<<'\t'<< v.E_interface << '\t'<<v.E_plus<<'\t'<<v.E_minus
		      <<'\t'<<v.plus_nodes<<'\t'<<v.minus_nodes<<'\n' ;
	ofile.close();
}

void CreateParams()
{
	string catalogue = "";
	string name = "params_";

	for (int i = 0; i < 10; i++) {
		ofstream ofile(catalogue + name + to_string(i));
		ofile << "model=qvoter_same" << '\n';
		ofile << "network=er" << '\n';
		ofile << "N=1000" << '\n';
		ofile << "c=0.5" << '\n';
		ofile << "k=4" << '\n';
		ofile << "q=1" << '\n';
		ofile << "t_max=20000000" << '\n';
		ofile << "p=0.0" << '\n';
		ofile << "folder=/home/joanna/workspace-cdt/network_sim/" << '\n';
		ofile << "filename=test_" + to_string(i) << '\n';
		ofile.close();
	}
}

// create files with params
void GenerateParams()
{
	string catalogue = "";
	string name = "params-";
	int Nvec[] = { 100,200,400,800,1600,3200 };
	int N;
	for (int j = 0; j < 6; j++) {
		N = Nvec[j];
		for (int p = 0; p <= 20; p++) {
			for (int i = 0; i < 1000; i++) {
				string ofile_name = catalogue + name + to_string(i) + "-" + to_string(p) + "-" + to_string(N);
				std::ofstream ofile(ofile_name);
				ofile << "model=qvoter_same" << '\n';
				ofile << "network=er" << '\n';
				ofile << "N=" + to_string(N) << '\n';
				ofile << "c=0.5" << '\n';
				ofile << "k=8" << '\n';
				ofile << "q=1" << '\n';
				ofile << "t_max=20000000" << '\n';
				ofile << "p=" + to_string((double(p)) / 20.0) << '\n';
				ofile << "folder=/home/joanna/workspace-cdt/network_sim/wyniki/same/" << '\n';
				ofile << "filename=" + to_string(i) << '\n';
				ofile.close();
			}
		}
	}
}
