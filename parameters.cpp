//Все, что связано с генерацией нового, случайного вектора
//1. выбор случайной координаты x для изменения
//2. генерация шага по выбранной координате. отклонения распределены
//по гауссу вокруг x с дисперсией dx
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "titles.h"

using namespace std;

Ran::Ran(unsigned long long int j) : v(4101842887655102017LL), w(1) {
	u = j^v; int64();
	v = u; int64();
	w = v; int64();
}
unsigned long long int Ran::int64() {
	u = u* 2862933555777941757LL + 7046029254386353087LL;
	v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
	w = 4294957665U* (w & 0xffffffff) + (w >> 32);
	unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
	return (x + v) ^ w;
}
double Ran::doub() {return 5.42101086242752217E-20* int64();}
unsigned int Ran::int32() {return (unsigned int) int64();}



Normaldev::Normaldev(double mmu, double ssig, unsigned long long i)
	: Ran(i), mu(mmu), sig(ssig) {}
double Normaldev::dev() {
	double u, v, x, y, q;
	do {
		u = doub();
		v = 1.7156* (doub() - 0.5);
		x = u - 0.449871;
		y = fabs(v) + 0.386595;
		q = pow(x, 2) + y* (0.19600*y - 0.25472*x);
	} while (q > 0.27597 && (q > 0.27846 || pow(v, 2) > -4. *log(u)*pow(u, 2)));
	return mu + sig*v/u;
}



Parameters::Parameters(int ranseed) : gau(0, 1, ranseed), data(){
	ifstream fin1("in/Parameters.txt"), fin2("in/Order.txt");
	ofstream fout("out/mcmc.txt");
	double q;
	string a, b;
	fin1 >> a;
	fin2 >> b;
	fin2 >> b;
	while(!fin1.eof()){
		if (a == b){ 
			fin1 >> q;
			x.push_back(q);
			propx.push_back(q);
			
			fin1 >> q;
			dx.push_back(q);

			fin1 >> q;
			xmin.push_back(q);

			fin1 >> q;
			xmax.push_back(q);
		}
		else cout << "Bad order in Parameters" << endl;
		fin1 >> a;
		fin2 >> b;
		fin2 >> b;
	}
	fin1.close();
	fin2.close();
	fout.close();
	FullSpectrum(x);
	data.update();
	data.chsq(chi2, prob);
}
void Parameters::propose(){
	int point, flag = 0;
	double delta;
	for(int i = 0; i < x.size(); i++) propx[i] = x[i];
	while(flag != 1){
		point = int(gau.doub()* x.size()); 
		delta = dx[point]* gau.dev();
		if (((propx[point] + delta > xmin[point]) && (propx[point] + delta < xmax[point])) && (dx[point] != 0.)) {
			flag = 1;
			propx[point] += delta;
		}
	}
	FullSpectrum(propx);
	data.update();
	data.chsq(pchi2, pprob);

}
void Parameters::update(){
	chi2 = pchi2;
	prob = pprob;
	for(int i = 0; i < x.size(); i++) x[i] = propx[i];
}
void Parameters::print(){
	ofstream fout("out/mcmc.txt", iostream::app);
	fout << 1 << "  " << prob << "  ";
	for(int i = 0; i < x.size(); i++) fout << x[i] << "  ";
	fout << endl;
	fout.close();
}

