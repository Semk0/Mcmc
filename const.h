#include <cmath>

const double Pi = atan(1)*4;					//пи 
const double h = 6.63*pow(10, -34);				//постоянная Планка h [Дж*с]
const double kb = 1.38*pow(10, -23);			//постояння Больцмана k [Дж/К]
const double c = 3*pow(10, 8);					//скорость света с [м/c]
const double secinyear = 31536000;				//секунд в году
const double Mpc = 3.085*pow(10, 22);			//метров в мегапарсеке
const double ae = 1.495*pow(10, 11);			//астрономическая единица в метрах
const double m_hyd = 1.66*pow(10, -27);			//масса атома водорода в кг
const int nuni = 1000;							//число делений радиуса
const int nm = 101;								//количество делений распределения звезд по массе
const int nz = 10000;							//количество делений z
const int ndust = 81;							//количество делений распределения пыли по размеру
const int nrad = 50;							//количество делений радиуса облака
const int nspec = 241;							
double	ksh;	
double	k1;
double	k3;
double	k4;
double	dm;
double	dr_cloud;
double	da; 
double runi;
double dr;


