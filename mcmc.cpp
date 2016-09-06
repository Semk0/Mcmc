#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include "titles.h"

using namespace std;


void mcmcstep(Parameters &param, double &rand) {
	double accept;
	accept = min(1., param.pprob/param.prob);
	if(rand < accept) {
		param.update();
		param.print();
	}
	else param.print();
}


int main(){
	double rand;
	time_t sec = time(NULL);
	Parameters param(sec);
	param.print();
	cout << param.chi2 << "   " << param.prob << endl;

	for(int i = 0; i < 10; i++) {
		param.propose();
		rand = param.gau.doub();
		mcmcstep(param, rand);
		cout << param.chi2 << "   " << param.prob << endl;
	}
	return 0;
}
