#include <vector>
#include <cmath>
using namespace std;



//spec
struct Spectrum;
struct Supp;
double zt(double, vector<double> &);
double tz(double, vector<double> &);
double rz(double, vector<double> &, Supp &);
double zr(double, vector<double> &, Supp &);
void Initialization(vector<double> &, Spectrum &, Supp &);
double LayerTemp(double, double, double, int, int, vector<double> &, Spectrum &, Supp &);
void StarSpectrum(double, double, double, vector<double> &, Spectrum &, Supp &);
void StarDustSpectrum(double, double, double, vector<double> &, Spectrum &, Supp &);
double StarCount(double, int, int, vector<double> &, Supp &);
void GalaxySpectrum(double, vector<double> &, Spectrum &, Supp &);
void Print(Spectrum &);
void FullSpectrum(vector<double> &);






//chisquare
double gammln(const double);
struct Gauleg18 {
	static const int ngau;
	static const double y[18];
	static const double w[18];
};
struct Gamma : Gauleg18{
	int ASWITCH;
	double EPS, FPMIN;
	double gln;
	Gamma();
	double gammp(const double, const double);
	double gammq(const double, const double); 
	double gser(const double, const double);
	double gcf(const double, const double); 
	double gammpapprox(double, double, int);
};
struct Datastr{
	int knstrn;
	vector<double> expx, expy, sigma, theor;
	Datastr();
	void update();
	void chsq(double &, double &);
};






//proposal
struct Ran {
	unsigned long long int u, v, w;
	Ran(unsigned long long int);
	unsigned long long int int64();
	double doub();
	unsigned int int32();
};
struct Normaldev : Ran {
	double mu, sig;
	Normaldev(double, double, unsigned long long int);
	double dev();
};
struct Parameters{
	double chi2, prob, pchi2, pprob; 
	vector<double> x, dx, xmin, xmax, propx;
	Normaldev gau;
	Datastr data;
	Parameters(int);
	void propose();
	void update();
	void print();
};





//mcmc
void mcmcstep(Parameters &, double &);
