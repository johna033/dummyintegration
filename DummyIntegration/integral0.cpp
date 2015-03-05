#include "stdafx.h"
#include "integral0.h"
#include <math.h>

double i0_fx(double x){
	return 4.5*cos(7*x)*exp((-2*x)/3) + 1.4*sin(1.5*x)*exp(-x/3) + 3.0;
}

double i0_px(double x){

	return pow((x - 2.1), -0.4);
}

double i0_mu0(double x){

	return (5.0/3.0)*pow((x - 2.1), 0.6);
}
double i0_mu1(double x){

	return (5.0/16.0)*pow((x - 2.1), 0.6)*(2*x + 7);
}
double i0_mu2(double x){

	return 5.0/208*pow((x - 21.0/10), 0.6)*(16*x*x + 42*x + 147);
}
double i0_mu3(double x){

	return (1.9036539387158786/7488)*pow((5*x - 21.0/2), 0.6)*(416*x*x*x + 1008*x*x + 2646*x + 9261);
}
double i0_mu4(double x){

	return (1.9036539387158786/28704)*pow((5*x - 21.0/2), 0.6)*(1248*x*x*x*x + 2912*x*x*x + 7056*x*x + 18522*x + 64827);
}
double i0_mu5(double x){

	return (1.9036539387158786/535808)*pow((5*x - 21.0/2), 0.6)*(19136*x*x*x*x*x + 43680*x*x*x*x+101920*x*x*x+246960*x*x+648270*x+2268945);
}