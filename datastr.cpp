//все, что связано с вычислением хи-квадрат 

using namespace std;
#include <cmath>
#include <iostream>
#include <fstream>
#include "spline.h"
#include "titles.h"

double gammln(const double xx) {
	int j;
	double x, tmp, y, ser;
	static const double cof[14] = {57.1562356658629235,
		-59.5979603554754912, 14.1360979747417471,
		-0.491913816097620199, .339946499848118887e-4,
		.465236289270485756e-4, -.983744753048795646e-4,
		.158088703224912494e-3, -.210264441724104883e-3,
		.217439618115212643e-3, -.164318106536763890e-3,
		.844182239838527433e-4, -.261908384015814087e-4,
		.368991826595316234e-5};
	if (xx <= 0) {
		cout << "bad arg in gammln" << endl;
		cin.get();}
	y = x = xx;
	tmp = x + 5.2421875;
	tmp = (x + 0.5)* log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0; j < 14; j++) ser += cof[j]/++y;
	return tmp + log(2.5066282746310005* ser/x);
}

const int Gauleg18::ngau = 18;
const double Gauleg18::y[18] = {0.0021695375159141994, 0.011413521097787704,
	0.027972308950302116, 0.051727015600492421, 0.082502225484340941,
	0.12007019910960293, 0.16415283300752470, 0.21442376986779355,
	0.27051082840644336, 0.33199876341447887, 0.39843234186401943,
	0.46931971407375483, 0.54413605556657973, 0.62232745288031077,
	0.70331500465597174, 0.78649910768313447, 0.87126389619061517,
	0.95698180152629142};
const double Gauleg18::w[18] = {0.0055657196642445571, 0.012915947284065419,
	0.020181515297735382, 0.027298621498568734, 0.034213810770299537,
	0.040875750923643261, 0.047235083490265582, 0.053244713977759692,
	0.058860144245324798, 0.064039797355015485, 0.068745323835736408,
	0.072941885005653087, 0.076598410645870640, 0.079687828912071670,
	0.082187266704339706, 0.084078218979661945, 0.085346685739338721,
	0.085983275670394821};


	Gamma::Gamma(){
		ASWITCH = 100;
		EPS = 1E-15;
		FPMIN = 1E-30;
	}
	
	double Gamma::gammp(const double a, const double x) {
		if (x < 0.0 || a <= 0.0) cout << "bad args in gammp";
		if (x == 0.0) return 0.0;
		else if (int(a) >= ASWITCH) return gammpapprox(a, x, 1);
		else if (x < a + 1.0) return gser(a, x);
		else return 1.0 - gcf(a, x);
	}
	double Gamma::gammq(const double a, const double x) {
		if (x < 0.0 || a <= 0.0) cout << "bad args in gammp";
		if (x == 0.0) return 1.0;
		else if (int(a) >= ASWITCH) return gammpapprox(a, x, 0);
		else if (x < a + 1.0) return 1.0 - gser(a, x);
		else return gcf(a, x);
	}
	double Gamma::gser(const double a, const double x) {
		double sum, del, ap;
		gln = gammln(a);
		ap = a;
		del = sum = 1.0/a;
		for(;;) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				return sum* exp(-x + a*log(x) - gln);
			}
		}
	}
	double Gamma::gcf(const double a, const double x) {
		int i;
		double an, b, c, d, del, h;
		gln = gammln(a);
		b = x + 1.0 - a;
		c = 1.0/FPMIN;
		d = 1.0/b;
		h = d;
		for (i = 1;;i++) {
			an = -i*(i - a);
			b += 2.0;
			d = an*d + b;
			if (fabs(d) < FPMIN) d = FPMIN;
			c= b + an/c;
			if (fabs(c) < FPMIN) c = FPMIN;
			d = 1.0/d;
			del = d*c;
			h *= del;
			if (fabs(del - 1.0) <= EPS) break;
		}
		return exp(-x + a*log(x) - gln)*h;
	}
	double Gamma::gammpapprox(double a, double x, int psig) {
		int j;
		double xu, t, sum, ans;
		double a1 = a - 1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
		gln = gammln(a);
		if (x > a1) xu = max(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
		else xu = max(0., min(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
		sum = 0;
		for (j = 0; j < ngau; j++) {
			t = x + (xu - x)*y[j];
			sum +=w[j]*exp(-(t - a1) + a1*(log(t) - lna1));
		}
		ans = sum*(xu - x)*exp(a1*(lna1 - 1.) - gln);
		return (psig?(ans > 0.0? 1.0 - ans:-ans):(ans >= 0.0? ans:1.0 + ans));
	}

Datastr::Datastr() {
	knstrn = 1;     //number of constraints
	ifstream fin("in/lower_limits.txt");
	double xx;
	while(!fin.eof()){
		fin >> xx;
		expx.push_back(xx);
		fin >> xx;
		expy.push_back(xx);
		fin >> xx;
		sigma.push_back(xx);		
	}
	fin.close();
}
void Datastr::update(){
	ifstream fin("in/Spectrum.txt");
	double xt;
	tk::spline spln;
	vector<double> xx, yy;
	theor.clear();
	fin >> xt;
	while(!fin.eof()){
		xx.push_back(xt);
		fin >> xt;
		yy.push_back(xt);
		fin >> xt;	
	}		
	fin.close();
	spln.set_points(xx, yy);
	for(int i = 0; i < expx.size(); i++) {
		xt = spln(expx[i]);
		theor.push_back(xt);
	}
}
void Datastr::chsq(double &chi2, double &prob){
	Gamma gam;
	int nbins = theor.size();
	double temp, df = nbins - knstrn;
	chi2 = 0.0;
	for(int j = 0; j < nbins; j++) {
		temp = theor[j] - expy[j];
		chi2 += temp*temp/sigma[j];
		}
	prob = gam.gammq(0.5*df, 0.5*chi2);
}

