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
#include <chrono>
#include <set>

using namespace std;

struct single_measure {
	int E_interface;
	int E_plus;
	int E_minus;
	int plus_nodes;
	int minus_nodes;
	int time;
};

void CreateGraphER(vector < set<int> > &network, int N, double k);
int CreateNodesState(int N, double c, vector<int> &nodes_states);
single_measure Measure(vector<int> &nodes_states, int N, vector < vector <int> > &network, int time);
//bool isConnected(vector < vector <int> > &network, int node1, int node2);
void RewireToRandom(vector < set<int> > &network, int basic_node, int old_neighbor, int N);
void RewireToSame(vector < set<int> > &network, int basic_node, int old_neighbor, 
	vector<int> &nodes_states, int inTheSameState);
bool CheckQpanel(vector<int> &qpanel, vector<int> &nodes_states, int q);
void DynamicsQVoterRandom(list<single_measure> &measure_list, int t_max,
	vector<set<int> > &network, vector<int> &nodes_states, int N, int q,
	double p);
void DynamicsQVoterSame(list<single_measure> &measure_list, int t_max,
	vector<set<int> > &network, vector<int> &nodes_states, int N, int q,
	double p);
/*void DynamicsQVoterSameAK(list<single_measure> &measure_list, int t_max,
		vector<set<int> > &network, vector<int> nodes_states, int N, int q,
		double p);*/
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

	auto start = std::chrono::system_clock::now();
#ifdef _DEBUG
	model = "qvoter_same";
	// ?? network = argv[2];
	N = 400;
	c = 0.5;
	k = 8;
	q = 1;
	p = 0.55;
	t_max = 10000;
	folder = "./";
	filename = "1";
#else

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
	//parameters in the order: model network N c k q p t_max folder filename
	model = argv[1];
	// ?? network = argv[2];
	N = atoi(argv[3]);
	c = atof(argv[4]);
	k = atof(argv[5]);
	q = atoi(argv[6]);
	p = atof(argv[7]);
	t_max = atoi(argv[8]);
	folder = argv[9];
	filename = argv[10];
#endif

	seed = (unsigned int) (time(NULL));
	if (argc > 11) {
		istringstream(argv[11]) >> seed;
		cout << "seed " << seed << endl;
	} else
		cout << "seed was chosen from the current time\n";
	srand(seed);

	cout << "seed " << seed << endl;
	cout << "params: model-" << model << "  N-" << N << "  c-" << c << "  k-"
			<< k << "  q-" << q << "  p-" << p << "  t_max-" << t_max
			<< "  filename-" << filename << "  folder-" << folder << "\n";

	list<single_measure> measure_list = { };

	outname = folder + model + "_q-" + to_string(q) + "_k-" + to_string(k)
			+ "_c-" + to_string(c) + "_N-" + to_string(N) + "_p-" + to_string(p)
			+ "_" + filename;
	cout << "output file name is " << outname << endl;

	vector<int> nodes_states(N);
	vector<set<int> > network(N, set<int>());

	CreateGraphER(network, N, k);
	cout << "CreateGraphER() done \n";
	CreateNodesState(N, c, nodes_states);
	cout << "CreateNodesState() done \n";

	if (model == "qvoter_same") {
		DynamicsQVoterSame(measure_list, t_max, network, nodes_states, N, q, p);
	} else if (model == "qvoter_random") {
		DynamicsQVoterRandom(measure_list, t_max, network, nodes_states, N, q,
				p);
	} else if (model == "qvoter_same_ak") {
		//DynamicsQVoterSameAK(measure_list,  t_max,network, nodes_states, N,  q,
		//		p);
	} else {
		cerr << "Wrong model name. Possible models: qvoter_same, qvoter_random, qvoter_same_ak"
				<< "\n";
		exit(-1);
	}

	SaveMeasure(measure_list, outname);

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";

	return 0;
}

// create a random Erdos Renyi graph with N nodes and mean degree k
void CreateGraphER(vector < set<int> > &network, int N, double k)
{
	double p=k/((double)N-1);
	double r;
	for (int i=0; i<N;i++){
		for(int j=i; j<N;j++){
			r =((double)rand()) / RAND_MAX;
			if(r<p){
				network[i].insert(j);
				network[j].insert(i);
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
single_measure Measure(vector<int> &nodes_states, int N, vector < set<int> > &network, int time)
{ //przemyslec co z tym wskaznikiem, int* E_int){
	int plus_nodes=0;
	int minus_nodes=0;
	int E_int=0;
	int E_plus=0;
	int E_minus=0;
	for (int i=0; i<N; i++){
		//size_t neighbor_number=network[i].size();
		if(nodes_states[i] == 1){
			plus_nodes++;
			//for(int j=0;j<neighbor_number;j++){
			for (auto x : network[i]) {
				//if(nodes_states[network[i][j]]==1)
				if (nodes_states[x]==1)
					E_plus++;
				else
					E_int++;
			}
		}else{
			minus_nodes++;
			//for(int j=0;j<neighbor_number;j++){
			for (auto x : network[i]) {
			//if(nodes_states[network[i][j]]==-1)
				if (nodes_states[x] == -1)
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
/*bool isConnected(vector < vector <int> > &network,int node1, int node2)
{
	bool connected=false;
	if (find(network[node1].begin(), network[node1].end(),node2)!=network[node1].end())
		connected=true;
	return connected;
}*/

//rewire basic_node link from old_neighbor to randomly chosen node in the network
void RewireToRandom(vector < set<int> > &network, int basic_node,int old_neighbor, int N )
{
	//TODO lack of ensuracne that there is no more not connected nodes in the network
	int new_neighbor=rand() % N;
	//while(isConnected(network, basic_node, new_neighbor))
	while (network[basic_node].count(new_neighbor)) // while new_neighbor is in the neighborhood 
		new_neighbor=rand() % N;
	//replace (network[basic_node].begin(), network[basic_node].end(),old_neighbor , new_neighbor);
	network[basic_node].erase(old_neighbor);
	network[basic_node].insert(new_neighbor);
	
	//vector<int>::iterator position = find(network[old_neighbor].begin(), network[old_neighbor].end(), basic_node);
	//network[old_neighbor].erase(position);
	network[old_neighbor].erase(basic_node);
	//network[new_neighbor].push_back(basic_node);
	network[new_neighbor].insert(basic_node);
}

// rewire basic_node link from old_neighbor to randomly chosen node in the network with the same state as basic_node
void RewireToSame(vector < set<int> > &network, int basic_node,int old_neighbor, vector<int> &nodes_states, int inTheSameState)
{
//very very slow probably
	int neighbor_in_the_same_state=0;
	int state_basic_nodes=nodes_states[basic_node];
	/*for(unsigned i=0;i<network[basic_node].size();i++){
		if(nodes_states[network[basic_node][i]]==state_basic_nodes)
			neighbor_in_the_same_state++;
	}*/
	for (auto x : network[basic_node]) {
		if (nodes_states[x] == state_basic_nodes)
			neighbor_in_the_same_state++;
	}
	int number_of_possible_nodes=inTheSameState-neighbor_in_the_same_state;
	if(number_of_possible_nodes>0){
		int rand_number=rand() % (number_of_possible_nodes);
		int j=0, new_neighbor=-1;
		while(j<=rand_number){
			new_neighbor++;
			//if(nodes_states[new_neighbor]==state_basic_nodes && !isConnected(network, basic_node, new_neighbor))
			if ((nodes_states[new_neighbor] == state_basic_nodes) && (network[basic_node].count(new_neighbor)==0))
				j++;
		}

		//replace (network[basic_node].begin(), network[basic_node].end(), old_neighbor, new_neighbor);
		network[basic_node].erase(old_neighbor);
		network[basic_node].insert(new_neighbor);
		//vector<int>::iterator position = find(network[old_neighbor].begin(), network[old_neighbor].end(), basic_node);
		//network[old_neighbor].erase(position);
		network[old_neighbor].erase(basic_node);
		//network[new_neighbor].push_back(basic_node);
		network[new_neighbor].insert(basic_node);
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
		vector<set<int> > &network, vector<int> &nodes_states, int N, int q,
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
			//size = network[basic_node].size();
			vector<int> v;
			v.assign(network[basic_node].begin(), network[basic_node].end());
			size = v.size();
			if (size > 0) {
				r1 = rand() % size;
				//qpanel[0] = network[basic_node][r1];
				qpanel[0] = v[r1];
				if (nodes_states[basic_node] != nodes_states[qpanel[0]]) {
					for (int i = 1; i < q; i++) {
						r1 = rand() % size;
						//qpanel[i] = network[basic_node][r1];
						qpanel[i] = v[r1];
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
		vector<set<int> > &network, vector<int> &nodes_states, int N, int q,
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
			//size = network[basic_node].size();
			vector<int> v;
			v.assign(network[basic_node].begin(), network[basic_node].end());
			size = v.size();
			if (size > 0) {
				r1 = rand() % size;
				//qpanel[0] = network[basic_node][r1];
				qpanel[0] = v[r1];
				if (nodes_states[basic_node] != nodes_states[qpanel[0]]) {
					for (int i = 1; i < q; i++) {
						r1 = rand() % size;
						//qpanel[i] = network[basic_node][r1];
						//
						qpanel[i] = v[r1];
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

/*
void DynamicsQVoterSameAK(list<single_measure> &measure_list, int t_max,
		vector<set<int> > &network, vector<int> nodes_states, int N, int q,
		double p){
cout<<"ak"<<'\n';
	int t = 0;
	vector<int> qpanel(q);
	int basic_node, r1, old_neighbor, inTheSameState,state,a;
	int number_of_plus_nodes = 0, number_of_minus_nodes = 0;
	for (int i = 0; i < N; i++) {
		if (nodes_states[i] == 1) {
			number_of_plus_nodes++;
		} else {
			number_of_minus_nodes++;
		}
	}
	size_t size;
	int neighbors_it, neighbor_counter;
	measure_list.push_back(Measure(nodes_states, N, network, t));

	while ((t < t_max) && (measure_list.back().E_interface > 0)) { //E_int>0){
		for (int microstep = 0; microstep < N; microstep++) {
			basic_node = rand() % N;
			size = network[basic_node].size();
			state=nodes_states[basic_node];
			a = 0;
			if (size > 0) {
				//for (int it = 0; it < size; it++) {
				for (auto x : network[basic_node]) {
					//if (nodes_states[network[basic_node][it]] != state)
					if (nodes_states[x] != state)
						a++;
				}
				if ((((double) rand()) / RAND_MAX)
						< (pow(((double) a) / ((double) size), q))) {
					r1 = rand() % a;
					neighbors_it = -1;
					neighbor_counter=-1;
					while (neighbor_counter < r1) {
						neighbors_it++;
						if (state != nodes_states[network[basic_node][neighbors_it]]) {
							neighbor_counter++;
						}
					}
					old_neighbor = network[basic_node][neighbors_it];
					if ((((double) rand()) / RAND_MAX) < p) {

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
						nodes_states[basic_node] = nodes_states[old_neighbor];
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
*/

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