#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "spline.h"
#include "const.h"
#include "titles.h"
using namespace std;

struct Spectrum{
	Spectrum();
	double l[nspec];
	double s[nspec];
	double g[nspec];
	double f[nspec];};
struct Supp{
	Supp(vector<double> &);
	double r[nm];
	double time[nm];
	double T[nm];
	double cmf[nm];				
	double sfr[nz];				
	double rn[nz];				
	double a[ndust];			
	double gra[ndust][nspec];
	double sil[ndust][nspec];};


//z from t
inline double zt(double t, vector<double> &x) {return pow(x[29]/x[28], 0.333)*pow(sinh(1.5*pow(x[29], 0.5)*t/x[27]), -0.666) - 1;}

//t from z
inline double tz(double z, vector<double> &x) {return 2/(3*pow(x[29], 0.5))*x[27]*log(pow(x[29]/(x[28]*pow(1 + z, 3)), 0.5) + pow(x[29]/(x[28]*pow(1 + z, 3)) + 1, 0.5));} 

//r from z
inline double rz(double z, vector<double> &x, Supp &s) {return s.rn[int(z/x[0]*(nz - 1))];}

//z from r
double zr(double r, vector<double> &x, Supp &s){
	int i_min = 0;
	int i_max = nz - 1;
	while (i_max - i_min > 1){
		if (s.rn[(i_min + i_max)/2] > r) i_max = (i_min + i_max)/2;
		else i_min = (i_min + i_max)/2;
	}
	return i_min*1.0/nz*x[0];
}

Spectrum::Spectrum(){
	ifstream fin1("in/Gra_81.txt");
	double w;
	string str;
	for(int h = 0; h < 11; h++) fin1 >> str;	
	for(int j = nspec - 1; j > -1; j--){						
 		fin1 >> l[j];
		l[j] = l[j]*pow(10, -6);		
		fin1 >> w;
		fin1 >> w;
		fin1 >> w;
	}
	fin1.close();
}
Supp::Supp(vector<double> &x){
	ifstream fin1, fin2;
	tk::spline spl;
	double w = 0;
	string str;    

	ksh = exp(-pow(x[5], 2)/(2*x[6]));	
	k1 = x[16]* pow(x[9], x[13] - x[12]);
	k3 = x[16]* pow(x[10], x[13] - x[14]);
	k4 = k3* pow(x[11], x[14] - x[15]);
	dm = (x[4] - x[3])/(nm - 1);
	dr_cloud = (x[19] - x[18])/nrad;
	da = (x[24] - x[23])/ndust;	
	
	fin1.open("in/Gra_81.txt");
	fin2.open("in/Sil_81.txt");
	for(int i = 0; i < ndust; i++){
		fin1 >> a[i];
		fin2 >> w;
		a[i] = a[i]* pow(10, -6);	
		for(int h = 0; h < 10; h++){		
			fin1 >> str;}
		for(int h = 0; h < 8; h++){			
			fin2 >> str;}	
		for(int j = nspec - 1; j > -1; j--){						
 			fin1 >> w;
			fin1 >> gra[i][j];
			fin1 >> w;
			fin1 >> w;
			fin2 >> w;
			fin2 >> sil[i][j];
			fin2 >> w;
			fin2 >> w;
		}}
	fin1.close();
	fin2.close();
//cmf	
	double f = 0;
	for(int i = 0; i < nm; i++){
		if (x[3] + i*dm < 1) cmf[i] = (1/(x[3] + i*dm))*
			exp(-pow((log10(x[3] + i*dm) - x[5]), 2)/(2*x[6]));
		else cmf[i] = ksh*pow(x[3] + i*dm, -x[7]);
	}
	for(int i = 0; i < nm - 1; i++) f += (x[3] + (i + 0.5)*dm)*(cmf[i + 1] + cmf[i])/2*dm;
	for(int i = 0; i < nm; i++) cmf[i] = cmf[i]/f;
//sfr
	fin1.open("in/sfr.txt");
	vector<double> Xc, Yc;
	fin1 >> w;
	while(!fin1.eof()){
		Xc.insert(Xc.begin(), tz(w, x));
		fin1 >> w;
		Yc.insert(Yc.begin(), pow(10, w));
		fin1 >> w;}
	spl.set_points(Xc,Yc);
	fin1.close();
	for(int j = 0; j < nz; j++) sfr[j] = spl(tz(x[0], x) + (tz(0, x) - tz(x[0], x))/(nz - 1)*j);
//rn
	double dx;
	rn[0] = 0;
	dx = x[0]/nz;
	for(int j = 1; j < nz; j++) rn[j] = rn[j - 1] + 
		c*secinyear*x[27]/Mpc*dx/pow(x[28]*pow(1 + j*dx, 3) + x[29], 0.5);
//red giants
	for(int j = 0; j < nm; j++){
		if(100*pow(x[3] + j*dm, x[8]) < 100) r[j] = 100* x[31]* pow(x[3] + j*dm, x[8]);
		else r[j] = 80 * x[31];
		time[j] = 0.09;
		T[j] = 4000;
	}
}


double LayerTemp(double T, double Rad, double r_layer, int j, int q, vector<double> &x, Spectrum &spec, Supp &s){
	double Iin = 0;			//падающая интенсивность
	double Iout = 0;		//исходящая интенсивность
	double dT = 1;			//точность решения
	double T_min = 1;		//начальная левая граница температуры
	double T_max = T;		//начальная правая граница температуры
	double C = 0;

	C = x[20]/((pow(x[23], x[26] + 1) - pow(x[24], x[26] + 1))/fabs(x[26] + 1) +
		(x[21]*x[33]*4*Pi)/(1.4*m_hyd*3)*(pow(x[24], x[26] + 4) - pow(x[23], x[26] + 4))/fabs(x[26] + 4));
	if(q == 0){
		for(int i = 0; i < nspec - 1; i++){
			Iin += 0.25*pow(Rad/r_layer, 2)* s.gra[j][i]
				*2*Pi*h*pow(c, 2)*pow(spec.l[i], -5)/(exp(h*c/(kb*T*spec.l[i])) - 1)*(spec.l[i + 1] - spec.l[i]);}
		while(T_max - T_min > dT){
			Iout = 0;
			for(int i = 0; i < nspec - 1; i++){
				Iout += s.gra[j][i]*2*Pi*h*pow(c, 2)*pow(spec.l[i], -5)/(exp(2*h*c/(kb*(T_max + T_min)*spec.l[i])) - 1)*(spec.l[i + 1] - spec.l[i]);}
			if (Iout < Iin) {
				T_min = (T_max + T_min)/2;}
			else {
				T_max = (T_max + T_min)/2;}}}
	else{
		for(int i = 0; i < nspec - 1; i++){
			Iin += 0.25*pow(Rad/r_layer, 2)* s.gra[j][i]
				*2*Pi*h*pow(c, 2)*pow(spec.l[i], -5)/(exp(h*c/(kb*T*spec.l[i])) - 1)*(spec.l[i + 1] - spec.l[i]);}
		while(T_max - T_min > dT){
			Iout = 0;
			for(int i = 0; i < nspec - 1; i++){
				Iout += s.sil[j][i]*2*Pi*h*pow(c, 2)*pow(spec.l[i], -5)/(exp(2*h*c/(kb*(T_max + T_min)*spec.l[i])) - 1)*(spec.l[i + 1] - spec.l[i]);}
			if (Iout < Iin) {
				T_min = (T_max + T_min)/2;}
			else {
				T_max = (T_max + T_min)/2;}}}
	return T_max;}

inline void StarSpectrum(double T, double r, double Rad, vector<double> &x, Spectrum &spec, Supp &s){
	for(int i = 0; i < nspec - 1; i++) spec.s[i] = pow(Rad/r, 2)* pow((1 + zr(r, x, s)), 2)* 
		2*Pi*h*pow(c, 2)*pow((spec.l[i] + spec.l[i + 1])/2, -5)/(exp(h*c*(1 + zr(r, x, s))/(kb*T*(spec.l[i] + spec.l[i + 1])/2)) - 1);
}
 
void StarDustSpectrum(double T, double r, double Rad, vector<double> &x, Spectrum &spec, Supp &s){
	double C = 0;						//константа в распределении пыли по размеру, определяется из условия нормировки на число частиц в единице объема
	double dS = 0;						//площадь поверхности частиц пыли в единичном объеме
	double Td;							//температура слойчика
	
	for(int i = 0; i < nspec; i++) {spec.s[i] = 0;}
	C = x[20]/((pow(x[23], x[26] + 1) - pow(x[24], x[26] + 1))/fabs(x[26] + 1) +
		(x[21]*x[33]*4*Pi)/(1.4*m_hyd*3)*(pow(x[24], x[26] + 4) - pow(x[23], x[26] + 4))/fabs(x[26] + 4));
	
	for(int i = 0; i < nrad; i++){
		for(int j = 0; j < ndust - 1; j++){
			/*if((x[23] < gra[j].a) && (gra[j].a < x[24])){
				Td = LayerTemp(T, Rad, (x[18] + i*dr_cloud)/Mpc, j, 0, x, spec, s);
				dS = 4*Pi*C*(pow(gra[j].a, -1.5) + pow(gra[j + 1].a, -1.5))/2*(gra[j + 1].a - gra[j].a)* 
					4/3*Pi*(pow(x[18] + (i + 1)*dr_cloud, 3) - pow(x[18] + i*dr_cloud, 3));
				dS = dS/pow(Mpc, 2); 		//перевод площади из метров в мегапарсеки
				for(int k = 0; k < nspec; k++){
					sds[k] += dS/(4*Pi*pow(r, 2))*gra[j].fabs[k]* pow((1+z[w]), 2)*
						2*Pi*h*pow(c, 2)*pow(gra[j].l[k], -5)/(exp(h*c*(1+z[w])/(kb*Td*gra[j].l[k])) - 1);}}}*/
			if((x[23] < s.a[j]) && (s.a[j] < x[24])){
				Td = LayerTemp(T, Rad, (x[18] + i*dr_cloud)/Mpc, j, 1, x, spec, s);
				dS = 4*Pi*C*(pow(s.a[j], x[26] + 2) + pow(s.a[j + 1], x[26] + 2))/2*(s.a[j + 1] - s.a[j])* 
					4/3*Pi*(pow(x[18] + (i + 1)*dr_cloud, 3) - pow(x[18] + i*dr_cloud, 3));
				dS = dS/pow(Mpc, 2); 		//перевод площади из метров в мегапарсеки
				for(int k = 0; k < nspec; k++){
					spec.s[k] += dS/(4*Pi*pow(r, 2))*s.sil[j][k]* pow((1+zr(r, x, s)), 2)*
						2*Pi*h*pow(c, 2)*pow(spec.l[k], -5)/(exp(h*c*(1 + zr(r, x, s))/(kb*Td*spec.l[k])) - 1);}}
		}}}

//число звезд данной массы
double StarCount(double r, int k, int type, vector<double> &x, Supp &s){
	double m = x[3] + (k + 0.5)*dm;
	int nc = 11;
	double dt;
	double ans = 0;
	double tfull, tini, tlife;
	switch (type){
	case 1:
		tlife = x[32]/pow(m, x[17]);
		if (tz(zr(r, x, s), x) - tz(x[0], x) > tlife){
			tfull = tlife;
			tini = tz(zr(r, x, s), x) - tlife;}
		else{
			tfull = tz(zr(r, x, s), x) - tz(x[0], x);
			tini = tz(x[0], x);}
		break;
	case 2:
		tlife = x[32]/pow(m, x[17]);
		if(tz(zr(r, x, s), x) - tz(x[0], x) > tlife){
			if(tz(zr(r, x, s), x) - tz(x[0], x) > tlife + (s.time[k] + s.time[k + 1])/2*tlife){
				tfull = (s.time[k] + s.time[k + 1])/2*tlife;
				tini = tz(zr(r, x, s), x) - tlife - (s.time[k] + s.time[k])/2*tlife;}
			else{
				tfull = tz(zr(r, x, s), x) - tlife;
				tini = tz(x[0], x);}}
		else {return 0;}
		break;
	case 3:
		tlife = x[25];
		if(tz(zr(r, x, s), x) - tz(x[0], x) > x[25]){
			tfull = x[25];
			tini = tz(zr(r, x, s), x) - x[25];}
		else{
			tfull = tz(zr(r, x, s), x) - tz(x[0], x);
			tini = tz(x[0], x);}
		break;
	default:
		break;}
	
	dt = tfull/(nc - 1);
	for(int j = 0; j < nc - 1; j++){ans += s.sfr[int((tini + (j + 0.5)*dt - tz(x[0], x))/(tz(0, x) - tz(x[0], x))*(nz - 1))]* pow(dr/(1 + zt(tini + (j + 0.5)*dt, x)), 3)* dt;}
	ans = ans* s.cmf[int(k)];
	return ans;}

//спектр галактики на расстоянии r
void GalaxySpectrum(double r, vector<double> &x, Spectrum &spec, Supp &s){
	double T, Rad, m;
	double number;
	for(int i = 0; i < nspec; i++) {spec.g[i] = 0;}
	for (int i = 0; i < nm - 1; i++){
		m = x[3] + (i + 0.5)*dm;
		Rad = x[31]*pow(m, x[8]);
		if((m < x[9]) && (m < x[4])) {T = pow(k1, 0.25)*x[30]*pow(m, 0.25*(x[12] - 2*x[8]));}
		if((m > x[9]) && (m < x[10]) && (m < x[4])) {T = pow(x[16], 0.25)*x[30]*pow(m, 0.25*(x[13] - 2*x[8]));}
		if((m > x[10]) && (m < x[11]) && (m < x[4])) {T = pow(k3, 0.25)*x[30]*pow(m, 0.25*(x[14] - 2*x[8]));}
		if((m > x[11]) && (m < x[4])) {T = pow(k4, 0.25)*x[30]*pow(m, 0.25*(x[15] - 2*x[8]));}

		StarSpectrum(T, r, Rad, x, spec, s);
		number = StarCount(r, i, 1, x, s);
		for(int j = 0; j < nspec; j++) spec.g[j] += spec.s[j]*number*dm;
		StarSpectrum((s.T[i] + s.T[i + 1])/2, r, (s.r[i] + s.r[i + 1])/2, x, spec, s);
		number = StarCount(r, i, 2, x, s);
		for(int j = 0; j < nspec; j++) spec.g[j] += spec.s[j]*number*dm;
		// if(m < 5){
		// 	number = StarCount(r, i, 3, x, s);
		// 	StarDustSpectrum(T, r, Rad, x, spec, s);
		// 	for(int j = 0; j < nspec; j++) spec.g[j] += spec.s[j]*number*dm;
		// }
	}}

//Вывод
void Print(Spectrum &spec){
	ofstream fout("in/Spectrum.txt");
	int i;
	for(i = 0; i < nspec - 2; i++) fout << (spec.l[i] + spec.l[i + 1])/2* pow(10, 10) << "   " << spec.f[i]* pow(10, 9)* (spec.l[i] + spec.l[i + 1])/2 << endl;
	fout  << (spec.l[i] + spec.l[i + 1])/2* pow(10, 10) << "   " << spec.f[i]* pow(10, 9)* (spec.l[i] + spec.l[i + 1])/2;

	fout.close();
}

void FullSpectrum(vector<double> &x){
	Spectrum spec;
	Supp s(x);
	runi = rz(x[0], x, s);
	dr = runi/nuni;
	int imin = rz(x[1], x, s)/dr + 1;
	int imax = rz(x[2], x, s)/dr;
	for(int j = 0; j < nspec; j++) spec.f[j] = 0;
	for(int i = imin; i < imax; i++){
		//unsigned int start_time =  clock(); // начальное время
		GalaxySpectrum(i*dr, x, spec, s);
		//unsigned int end_time = clock(); 
		//cout << i << "    " << (end_time - start_time)/1000. << endl;
		for(int j = 0; j < nspec; j++) spec.f[j] += spec.g[j]*pow(i*dr, 2)*pow(dr, -2);
	}
	Print(spec);
	//cout << "FullSpec ok" << endl;
}



