#include <iostream>
#include <fstream>

using namespace std;


void create_params() {
	string catalogue = "";
	string name = "params-";
	int Nvec[]={100,200,400,800,1600,3200};
	int N;
	for(int j=0; j<6;j++){
		N=Nvec[j];
		for (int p=0; p<=20;p++){
			for (int i = 0; i < 1000; i++) {
				std::ofstream file(catalogue + name + to_string(i)+"-"+to_string(p)+"-"+to_string(N));
				file << "model=qvoter_same" << '\n';
				file << "network=er" << '\n';
				file << "N="+to_string(N) << '\n';
				file << "c=0.5" << '\n';
				file << "k=8" << '\n';
				file << "q=1" << '\n';
				file << "t_max=20000000" << '\n';
				file << "p="+to_string((double(p))/20.0) << '\n';
				file << "folder=/home/joanna/workspace-cdt/network_sim/wyniki/same/" << '\n';
				file << "filename=" + to_string(i) << '\n';
			}
		}
	}
}

int main(int argc, char **argv)
{
	create_params();	

	return 0;
}
