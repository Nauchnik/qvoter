#include "network_sim.h"

#include <sstream>

#ifdef _WIN32
#include "dirent.h"
#else
#include <dirent.h>
#endif

network_simiulation_sequential::network_simiulation_sequential() :
	t0(0),
	t_max(1000000000),
	folder(""),
	filename(""),
	seed(-1),
	realization(0),
	status(-1),
	start_time(0.0),
	previous_launches_count(0)
{}

void network_simiulation_sequential::ReadParams(const int argc, char **argv)
{
	//parameters in the order: model network_type N c k q p t_max folder filename
	model = argv[1];
	network_type = argv[2];
	N = atoi(argv[3]);
	c = atof(argv[4]);
	k = atof(argv[5]);
	q = atoi(argv[6]);
	p = atof(argv[7]);
	t0 = atoi(argv[8]);
	t_max = atoi(argv[9]);
	folder = argv[10];
	filename = argv[11];

	if (argc > 12) {
		istringstream(argv[12]) >> seed;
		cout << "seed " << seed << endl;
	}
	else {
		seed = (unsigned int)(time(NULL));
		cout << "seed was chosen from the current time\n";
	}

	cout << "seed " << seed << endl;
	cout << "params: model-" << model << "  N-" << N << "  c-" << c << "  k-"
		<< k << "  q-" << q << "  p-" << p << " t0-" << t0 << "  t_max-" << t_max
		<< "  filename-" << filename << "  folder-" << folder << "\n";
}

void network_simiulation_sequential::GetOutputName() 
{
	stringstream sstream;
	sstream << folder << model << "_q-" << q << "_k-" << k << "_c-" << c << "_N-" << N << "_p-" << p;
	if (filename != "")
		sstream << "_" << filename;
	sstream << "_r-" << realization << "_seed-" << seed;
	
	outname = sstream.str();
	if (verbosity > 0)
		cout << "output file name is " << outname << endl;
}

void network_simiulation_sequential::Init()
{
	//ReadSimulationState("qvoter_same_q-2_k-8_c-0.5_N-10_p-0_test_r--1387736016_seed-1_state_in_step_2");
	srand(seed);
	nodes_states.resize(N);
	network.resize(N);
}

void network_simiulation_sequential::LaunchSimulation()
{

	if (model == "qvoter_same") {
		DynamicsQVoterSame(t0, measure_list, t_max, network, nodes_states, N, q, p);
	}
	else if (model == "qvoter_random") {
		DynamicsQVoterRandom(t0, measure_list, t_max, network, nodes_states, N, q,
			p);
	}
	else if (model == "qvoter_same_ak") {
		DynamicsQVoterSameAK(t0, measure_list, t_max, network, nodes_states, N, q,
			p);
	}
	else {
		cerr << "Wrong model name. Possible models: qvoter_same, qvoter_random, qvoter_same_ak"
			<< "\n";
		exit(-1);
	}
}

// create a random Erdos Renyi graph with N nodes and mean degree k
void network_simiulation_sequential::CreateGraphER()
{
	double p = k / ((double)N - 1);
	double r;
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) { //i+1 - no self edges
			r = ((double)rand()) / RAND_MAX;
			if (r < p) {
				network[i].push_back(j);
				network[j].push_back(i);
			}
		}
	}
}

//adds states to nodes, c*N nodes with plus state, (1-c)*N with -1 state
int network_simiulation_sequential::CreateNodesState()
{
	int bound = (int)round(c*(double)N);
	for (int i = 0; i < bound; i++)
		nodes_states[i] = 1;
	for (int i = bound; i < N; i++)
		nodes_states[i] = -1;
	return 0;
}

//measure number of links in each state and number of nodes in each state
single_measure network_simiulation_sequential::Measure(vector<int> &nodes_states, int N, vector < vector <int> > &network, int time)
{ //przemyslec co z tym wskaznikiem, int* E_int){
	int plus_nodes = 0;
	int minus_nodes = 0;
	int E_int = 0;
	int E_plus = 0;
	int E_minus = 0;
	for (int i = 0; i < N; i++) {
		size_t neighbor_number = network[i].size();
		if (nodes_states[i] == 1) {
			plus_nodes++;
			for (int j = 0; j < neighbor_number; j++) {
				if (nodes_states[network[i][j]] == 1)
					E_plus++;
				else
					E_int++;
			}
		}
		else {
			minus_nodes++;
			for (int j = 0; j < neighbor_number; j++) {
				if (nodes_states[network[i][j]] == -1)
					E_minus++;
				else
					E_int++;
			}
		}
	}
	single_measure cokolwiek;
	cokolwiek.plus_nodes = plus_nodes;
	cokolwiek.minus_nodes = minus_nodes;
	cokolwiek.E_interface = E_int / 2;
	cokolwiek.E_plus = E_plus / 2;
	cokolwiek.E_minus = E_minus / 2;
	cokolwiek.time = time;

	return cokolwiek;
}

//check if two nodes are connected
bool network_simiulation_sequential::isConnected(vector < vector <int> > &network, int node1, int node2)
{
	bool connected = false;
	if (find(network[node1].begin(), network[node1].end(), node2) != network[node1].end())
		connected = true;
	return connected;
}

//rewire basic_node link from old_neighbor to randomly chosen node in the network
void network_simiulation_sequential::RewireToRandom(vector < vector <int> > &network, int basic_node, int old_neighbor, int N)
{
	//TODO lack of ensuracne that there is no more not connected nodes in the network
	int new_neighbor = rand() % N;
	while (isConnected(network, basic_node, new_neighbor))
		new_neighbor = rand() % N;
	replace(network[basic_node].begin(), network[basic_node].end(), old_neighbor, new_neighbor);
	vector<int>::iterator position = find(network[old_neighbor].begin(), network[old_neighbor].end(), basic_node);
	network[old_neighbor].erase(position);
	network[new_neighbor].push_back(basic_node);
}

// rewire basic_node link from old_neighbor to randomly chosen node in the network with the same state as basic_node
void network_simiulation_sequential::RewireToSame(vector < vector <int> > &network, int basic_node, int old_neighbor, vector<int> &nodes_states, int inTheSameState)
{
	//very very slow probably
	int neighbor_in_the_same_state = 0;
	int state_basic_nodes = nodes_states[basic_node];
	for (unsigned i = 0; i < network[basic_node].size(); i++) {
		if (nodes_states[network[basic_node][i]] == state_basic_nodes)
			neighbor_in_the_same_state++;
	}
	int number_of_possible_nodes = inTheSameState - neighbor_in_the_same_state-1;//minus 1 - we cannot add self edges
	if (number_of_possible_nodes > 0) {
		int rand_number = rand() % (number_of_possible_nodes);
		int j = 0, new_neighbor = -1;
		while (j <= rand_number) {
			new_neighbor++;
			if (nodes_states[new_neighbor] == state_basic_nodes && !isConnected(network, basic_node, new_neighbor) &&new_neighbor!=basic_node)
				j++;
		}

		replace(network[basic_node].begin(), network[basic_node].end(), old_neighbor, new_neighbor);
		vector<int>::iterator position = find(network[old_neighbor].begin(), network[old_neighbor].end(), basic_node);
		network[old_neighbor].erase(position);
		network[new_neighbor].push_back(basic_node);
	}
}

bool network_simiulation_sequential::CheckQpanel(vector<int> &qpanel, vector<int> &nodes_states, int q)
{
	bool identical = true;
	int state = nodes_states[qpanel[0]];
	for (int i = 1; i < q; i++) {
		if (nodes_states[qpanel[i]] != state) {
			identical = false;
			break;
		}
	}
	return identical;
}

//whole random qvoter model dynamics
void network_simiulation_sequential::DynamicsQVoterRandom(int t0, list<single_measure> &measure_list, int t_max,
	vector<vector<int> > &network, vector<int> &nodes_states, int N, int q,
	double p)
{
	int t = t0;
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
						if ((((double)rand()) / RAND_MAX) < p) {
							r1 = rand() % q;
							old_neighbor = qpanel[r1];

							RewireToRandom(network, basic_node, old_neighbor,
								N);
						}
						else {
							nodes_states[basic_node] = nodes_states[qpanel[0]];
						}
					}
				}
			}
		}
		t++;
		measure_list.push_back(Measure(nodes_states, N, network, t)); //, &E_int));
		if ( ((t % 100) == 0) && (verbosity > 0) )
			cout << t << "\n";
	}
}

//whole same qvoter model dynamics
void network_simiulation_sequential::DynamicsQVoterSame(int t0, list<single_measure> &measure_list, int t_max,
	vector<vector<int> > &network, vector<int> &nodes_states, int N, int q,
	double p)
{
	int t = t0;
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
						if ((((double)rand()) / RAND_MAX) < p) {
							r1 = rand() % q;
							old_neighbor = qpanel[r1];

							if (nodes_states[basic_node] == 1) {
								inTheSameState = number_of_plus_nodes;
							}
							else {
								inTheSameState = number_of_minus_nodes;
							}
							RewireToSame(network, basic_node, old_neighbor,
								nodes_states, inTheSameState);
						}
						else {
							if (nodes_states[basic_node] == 1) {
								number_of_plus_nodes--;
								number_of_minus_nodes++;
							}
							else {
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
		if ((verbosity > 0) && (t % 100 == 0))
			cout << t << "\n";
#ifdef _MPI
		if (MPI_Wtime() - start_time >= MAX_SOLVING_TIME_SEC) {
			status = -2; // interrupted
			SaveSimulationState(t, model, number_of_plus_nodes, number_of_minus_nodes);
			break;
		}
#endif
	}
	if (status != -2) // if not interrupted
		status = 1; // then solved
}

void network_simiulation_sequential::DynamicsQVoterSameAK(int t0, list<single_measure> &measure_list, int t_max,
	vector<vector<int> > &network, vector<int> nodes_states, int N, int q,
	double p) 
{
	if (verbosity > 0)
		cout << "ak" << '\n';
	int t = t0;
	vector<int> qpanel(q);
	int basic_node, r1, old_neighbor, inTheSameState, state, a;
	int number_of_plus_nodes = 0, number_of_minus_nodes = 0;
	for (int i = 0; i < N; i++) {
		if (nodes_states[i] == 1) {
			number_of_plus_nodes++;
		}
		else {
			number_of_minus_nodes++;
		}
	}
	int size, neighbors_it, neighbor_counter;
	measure_list.push_back(Measure(nodes_states, N, network, t));

	while ((t < t_max) && (measure_list.back().E_interface > 0)) { //E_int>0){
		for (int microstep = 0; microstep < N; microstep++) {
			basic_node = rand() % N;
			size = network[basic_node].size();
			state = nodes_states[basic_node];
			a = 0;
			if (size > 0) {
				for (int it = 0; it < size; it++) {
					if (nodes_states[network[basic_node][it]] != state) {
						a++;
					}
				}
				if ((((double)rand()) / RAND_MAX)
					< (pow(((double)a) / ((double)size), q))) {
					r1 = rand() % a;
					neighbors_it = -1;
					neighbor_counter = -1;
					while (neighbor_counter < r1) {
						neighbors_it++;
						if (state
							!= nodes_states[network[basic_node][neighbors_it]]) {
							neighbor_counter++;
						}
					}
					old_neighbor = network[basic_node][neighbors_it];
					if ((((double)rand()) / RAND_MAX) < p) {

						if (nodes_states[basic_node] == 1) {
							inTheSameState = number_of_plus_nodes;
						}
						else {
							inTheSameState = number_of_minus_nodes;
						}
						RewireToSame(network, basic_node, old_neighbor,
							nodes_states, inTheSameState);
					}
					else {
						if (nodes_states[basic_node] == 1) {
							number_of_plus_nodes--;
							number_of_minus_nodes++;
						}
						else {
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
		if ( ((t % 100) == 0) && (verbosity > 0) )
			cout << t << "\n";
	}

}

// write the measure list to the output file
void network_simiulation_sequential::SaveMeasure()
{
	string filename = outname;
	if (status == -2)
		filename += "_interrupted";
	ofstream ofile(filename);
	ofile << "#time\tE_interface\tE_plus\tE_minus\tplus_nodes\tminus_nodes\n";
	for (const auto& v : measure_list)
		ofile << v.time << '\t' << v.E_interface << '\t' << v.E_plus << '\t' << v.E_minus
		<< '\t' << v.plus_nodes << '\t' << v.minus_nodes << '\n';
	ofile.close();
}

void network_simiulation_sequential::CreateParams()
{
	string catalogue = "";
	string name = "params_";

	for (int i = 0; i < 10; i++) {
		ofstream ofile(catalogue + name + to_string((long long)i));
		ofile << "model=qvoter_same" << '\n';
		ofile << "network=er" << '\n';
		ofile << "N=1000" << '\n';
		ofile << "c=0.5" << '\n';
		ofile << "k=4" << '\n';
		ofile << "q=1" << '\n';
		ofile << "t_max=20000000" << '\n';
		ofile << "p=0.0" << '\n';
		ofile << "folder=/home/joanna/workspace-cdt/network_sim/" << '\n';
		ofile << "filename=test_" + to_string((long long)i) << '\n';
		ofile.close();
	}
}

// create files with params
void network_simiulation_sequential::GenerateParams()
{
	string catalogue = "";
	string name = "params-";
	int Nvec[] = { 100,200,400,800,1600,3200 };
	int N;
	for (int j = 0; j < 6; j++) {
		N = Nvec[j];
		for (int p = 0; p <= 20; p++) {
			for (int i = 0; i < 1000; i++) {
				string ofile_name = catalogue + name + to_string((long long)i) + "-" 
					+ to_string((long long)p) + "-" + to_string((long long)N);
				ofstream ofile(ofile_name);
				ofile << "model=qvoter_same" << '\n';
				ofile << "network=er" << '\n';
				ofile << "N=" + to_string((long long)N) << '\n';
				ofile << "c=0.5" << '\n';
				ofile << "k=8" << '\n';
				ofile << "q=1" << '\n';
				ofile << "t_max=20000000" << '\n';
				ofile << "p=" + to_string((long double)(p / 20.0)) << '\n';
				ofile << "folder=/home/joanna/workspace-cdt/network_sim/wyniki/same/" << '\n';
				ofile << "filename=" + to_string((long long)i) << '\n';
				ofile.close();
			}
		}
	}
}


void network_simiulation_sequential::SaveSimulationState(int step, string model, int number_of_plus_nodes, int number_of_minus_nodes)
{
	ofstream ofile(outname+"_state_in_step_"+to_string((long long)step));
	if(model=="qvoter_same" || model=="qvoter_same_ak"){
		ofile << "#model\tq\tp\tN\tnumber_of_plus_nodes\tnumber_of_minus_nodes\tstep"<< '\n';
		ofile << model<<"\t"<<q<<"\t"<<p<<"\t"<<N<<"\t"<<number_of_plus_nodes<<"\t"<<number_of_minus_nodes<<"\t"<<step<< '\n';
	}else{
		ofile << "#model\tq\tp\tN\tnumber_of_plus_nodes\tnumber_of_minus_nodes\tstep"<< '\n';
		ofile << model<<"\t"<<q<<"\t"<<p<<"\t"<<N<<"\t"<<"xx"<<"\t"<<"xx"<<"\t"<<step<< '\n';
	}
	ofile<<"#node"<<"\t"<<"state"<<"\n";
	for (int i=0; i< nodes_states.size();i++){
		ofile<<i<<"\t"<<nodes_states[i]<<"\n";
	}
	ofile<<"edges"<<"\n";
	if (verbosity > 0)
		cout<<network.size();
	for(int i=0;i<network.size();i++){
		for(int j=0; j<network[i].size();j++){
			ofile<<i<<"\t"<<network[i][j]<<"\n";
		}
	}
	ofile.close();
}

vector<string> network_simiulation_sequential::split(const string& s, char delimiter)
{
   vector<string> tokens;
   string token;
   istringstream tokenStream(s);
   while (getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

void network_simiulation_sequential::ReadSimulationState(string inname)
{
	int number_of_plus_nodes, number_of_minus_nodes;
	string line;
	ifstream infile(inname);

	if (!infile.is_open()) {
		cerr << "Unable to open the file: " << inname << "\n";
		exit(-1);
	}

	getline(infile, line);
	getline(infile, line);
	vector < string > arr = split(line, '\t');
	model = arr[0];
	q = stoi(arr[1]);
	p = stof(arr[2]);
	N = stoi(arr[3]);

	if (model == "qvoter_same" || model == "qvoter_same_ak") {
		number_of_plus_nodes = stoi(arr[4]);
		number_of_minus_nodes = stoi(arr[5]);
	}
	t0 = stoi(arr[6]);
	cout << "model: " << model << "\t" << "q: " << q << "\t" << "p: " << p
			<< "\t" << "N: " << N << "\n";
	getline(infile, line);
	nodes_states.resize(N);
	network.resize(N);

	for (int i = 0; i < N; i++) {
		getline(infile, line);
		arr = split(line, '\t');
		nodes_states[stoi(arr[0])] = stoi(arr[1]);
	}
//		for (int i = 0; i < N; i++) {
//			cout << i << "   " << nodes_states[i] << "\n";
//		}
	getline(infile, line);
	while (getline(infile, line)) {
		arr = split(line, '\t');
		network[stoi(arr[0])].push_back(stoi(arr[1]));
	}

//		for (int i = 0; i < network.size(); i++) {
//			for (int j = 0; j < network[i].size(); j++) {
//				cout << i << "   " << network[i][j] << "\n";
//			}
//		}

	infile.close();
}

// execute command via system process
string network_simiulation_sequential::exec(string cmd_str)
{
	string result = "";
	char* cmd = new char[cmd_str.size() + 1];
	for (unsigned i = 0; i < cmd_str.size(); i++)
		cmd[i] = cmd_str[i];
	cmd[cmd_str.size()] = '\0';
#ifdef _WIN32
	FILE* pipe = _popen(cmd, "r");
#else
	FILE* pipe = popen(cmd, "r");
#endif
	delete[] cmd;
	if (!pipe) return "ERROR";
	char buffer[128];
	while (!feof(pipe)) {
		if (fgets(buffer, 128, pipe) != NULL)
			result += buffer;
	}
#ifdef _WIN32
	_pclose(pipe);
#else
	pclose(pipe);
#endif
	return result;
}

void network_simiulation_sequential::getdir(string dir, vector<string> &files)
{
	DIR *dp;
	string cur_name;
	struct dirent *dirp;
	if ((dp = opendir(dir.c_str())) == NULL) {
		cerr << endl << "Error in opening folder " << dir << endl;
		exit(-1);
	}
	while ((dirp = readdir(dp)) != NULL) {
		cur_name = string(dirp->d_name);
		if (cur_name[0] != '.')
			files.push_back(cur_name);
	}
	closedir(dp);
}

string network_simiulation_sequential::FindStateFileName()
{
	string result = "";
	vector<string> files_names = vector<string>();
#ifdef WIN32
	string cur_path = "./";
#else
	string cur_path = exec("echo $PWD");
#endif
	cur_path.erase(remove(cur_path.begin(), cur_path.end(), '\r'), cur_path.end());
	cur_path.erase(remove(cur_path.begin(), cur_path.end(), '\n'), cur_path.end());
	getdir(cur_path, files_names);
	int max_state_in_step = -1;
	for (auto file_name : files_names) {
		size_t pos1 = file_name.find(outname);
		size_t pos2 = file_name.find("state_in_step_");
		if ((pos1 != string::npos) && (pos2 != string::npos)) {
			previous_launches_count++;
			string sstr = file_name.substr(pos2 + 14);
			int cur_state_in_step = -1;
			istringstream(sstr) >> cur_state_in_step;
			if ((cur_state_in_step > 0) && (cur_state_in_step > max_state_in_step)) {
				max_state_in_step = cur_state_in_step;
				result = file_name;
			}
		}
	}
	return result;
}

