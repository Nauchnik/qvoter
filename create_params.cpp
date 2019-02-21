#include <iostream>
#include <fstream>

using namespace std;


void create_params() {
	std::ofstream file("params.txt");
	string catalogue = "";
	string name = "";
	int Nvec[]={800};
	int N;
	for(int j=0; j<1;j++){
		N=Nvec[j];
		for (int p=0; p<=20;p++){

			for (int i = 1; i < 1000; i++) {
				file<<"./netsim"<<'\t';
				file << "qvoter_same_ak" << '\t';
				file << "er" << '\t';
				file << to_string(N) << '\t';
				file << "0.5" << '\t';
				file << 8 << '\t'; //k
				file << 1 << '\t'; //q
				file << to_string((double(p))/20.0) << '\t'; //p
				file << 2000000000 << '\t';
				file << "/home/joanna/workspace-cdt/network_sim/wyniki/same_ak/" << '\t';
				file << name<<i<< '\n';
			}
		}
	}
}

int main(int argc, char **argv)
{
	create_params();	

	return 0;
}
