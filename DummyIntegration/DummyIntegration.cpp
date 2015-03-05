// DummyIntegration.cpp: 
//

#include "stdafx.h"
#include "sillymatrix.h"
#include "integral0.h"
#include <math.h>
#include <iostream>
using namespace std;

#define CARDANO 0
#define NEWTON 1
#define NEWTON_EPSILON 1e-9

#define M_PI (3.141592653589793)
#define M_2PI (2.0*M_PI)

typedef enum _qType { Q_NEWTON, Q_GAUSS} Q_TYPE;

typedef double (*calc_func)(double);

class dummyIntegral{
private:
	calc_func fx;
	calc_func px;
	calc_func* mu;
	int musCount;
	dummyIntegral(){};
	dummyIntegral(dummyIntegral& A){};
	double* solveCubicNewton(double* coeffs, double initial);
	double* solveCubicCardano(double* coeffs);
	double* buildThreePointGaussianQuadratic(double a, double b);
	double* buildThreePointNewton(double a, double b, double* nodes);
	double* buildApproximationMesh(double a, double b, int numberOfNodes);
	double integrateWithComposedQuadraticNewton(double a, double b, int numberOfSteps);
	double integrateWithComposedQuadraticGauss(double a, double b, int numberOfSteps);
	double processAitken(double a, double b, Q_TYPE quadraticType);
public:
	dummyIntegral(calc_func fx, calc_func px, calc_func* mu, int musCount):fx(fx), px(px), mu(mu), musCount(musCount){};
	~dummyIntegral(){};
	double integrateWithComposedQuadratic(double a, double b, Q_TYPE quadraticType);
	double integrateWithComposedQuadratic(double a, double b, double epsilon, Q_TYPE quadraticType);
	double processAitken(double a, double b, double epsilon, Q_TYPE quadraticType);
	double integrateWithOptimalStepThroughAitken(double a, double b, double epsilon, Q_TYPE quadraticType);
	double integrateWithOptimalStep(double a, double b, double epsilon, Q_TYPE quadraticType);
};
double* dummyIntegral::solveCubicCardano(double* coeffs){
	double q = (coeffs[2]*coeffs[2] - 3*coeffs[1])/9.0;
	double r = (2.0*coeffs[2]*coeffs[2]*coeffs[2] - 9.0*coeffs[2]*coeffs[1]+27.0*coeffs[0])/54.0;
	if(r*r < q*q*q){
		cout << fabs(r*r-q*q*q) << endl;
		double* roots = new double[3];
		double t = acos(r/sqrt(q*q*q));
		double u = coeffs[2]/3.0;
		roots[0] = -2*sqrt(q)*cos(t/3.0) - u;
		roots[1] = -2*sqrt(q)*cos((t+M_PI)/3.0) - u;
		roots[2] = -2*sqrt(q)*cos((t-M_2PI)/3.0) - u;
		return roots;
	}
	return NULL;
}

double* dummyIntegral::buildApproximationMesh(double a, double b, int numberOfNodes){

	double* net = new double[numberOfNodes];
	double step = (b - a)/numberOfNodes;
	cout <<"Step size: " << step << endl;
	for(int i = 0; i < numberOfNodes; i++){
		net[i] = a + i*step;
	}
	net[numberOfNodes - 1] = b;
	return net;
}
double* dummyIntegral::solveCubicNewton(double* coeffs, double initial){
	double* roots = new double[3];
	double approx = initial;
		while(1){
			roots[0] = approx - (approx*approx*approx + coeffs[2]*approx*approx + coeffs[1]*approx + coeffs[0])/(3*approx*approx + 2*coeffs[2]*approx + coeffs[1]);
			if(fabs(roots[0] - approx) <= NEWTON_EPSILON)
				break;
			approx = roots[0];
		}
	double b,c;
	b = (coeffs[2] + roots[0]);
	c = (coeffs[1] + roots[0]*b);
	double quadratic = b*b - 4*c;
	//cout << coeffs[0]+c*roots[0] << endl;
	if(quadratic <= NEWTON_EPSILON){
		delete[] roots;
		return NULL;
	}
	roots[1] = (-b+sqrt(quadratic))/2.0;
	roots[2] = (-b-sqrt(quadratic))/2.0;
	return roots;

}
double* dummyIntegral::buildThreePointGaussianQuadratic(double a, double b){
	double* weightMomentums = new double[musCount];
	for(int i = 0; i < musCount; i++){
		weightMomentums[i] = mu[i](b) - mu[i](a);
	}
	SillyMatrix* matrix = new SillyMatrix(3, 3);
	(*matrix)(0, 0) = weightMomentums[0]; (*matrix)(0, 1) = weightMomentums[1]; (*matrix)(0, 2) = weightMomentums[2];
	(*matrix)(1, 0) = weightMomentums[1]; (*matrix)(1, 1) = weightMomentums[2]; (*matrix)(1, 2) = weightMomentums[3];
	(*matrix)(2, 0) = weightMomentums[2]; (*matrix)(2, 1) = weightMomentums[3]; (*matrix)(2, 2) = weightMomentums[4];
	double* b1 = new double[3];
	b1[0] = -1*weightMomentums[3];
	b1[1] = -1*weightMomentums[4];
	b1[2] = -1*weightMomentums[5];
	int opCount = 0;
	double* cubicCoeffs = matrix->solveSystem(matrix, b1, &opCount);
	if(cubicCoeffs == NULL){
		delete[] weightMomentums;
		delete matrix;
		delete[] b1;
		throw "Can't solve SLAE for cubic coeffs!";
		return NULL;
	}
	int type = NEWTON;
	//cout << "Polynom : x^3+" << cubicCoeffs[2] << "*x^2+" << cubicCoeffs[1] << "*x+" << cubicCoeffs[0] << " = 0" << endl;
	//cout << "Solve using Cardano(0) or Newton(1) :"; cin >> type;
	double* roots = NULL;
	if(type == NEWTON){
		roots = this->solveCubicNewton(cubicCoeffs, b);
		//cout << roots[0] << " " << roots[1] << " " << roots[2] << endl;
	}else{
		roots = this->solveCubicCardano(cubicCoeffs);
		//cout << roots[0] << " " << roots[1] << " " << roots[2] << endl;
		//cin.get();
	}
	delete[] cubicCoeffs;
	if(roots != NULL){
		if(roots[0] > b || roots[0] < a || roots[1] > b || roots[1] < a || roots[2] > b || roots[2] < a){	
		delete[] weightMomentums;
		delete matrix;
		delete[] b1; 
		delete[] roots;
		throw "Roots out of range!";
		return NULL;
		}
	(*matrix)(0, 0) = 1.0; (*matrix)(0, 1) = 1.0; (*matrix)(0, 2) = 1.0;
	(*matrix)(1, 0) = roots[0]; (*matrix)(1, 1) = roots[1]; (*matrix)(1, 2) = roots[2];
	(*matrix)(2, 0) = roots[0]*roots[0]; (*matrix)(2, 1) = roots[1]*roots[1]; (*matrix)(2, 2) = roots[2]*roots[2];
	b1[0] = weightMomentums[0];
	b1[1] = weightMomentums[1];
	b1[2] = weightMomentums[2];
	double* bigA = matrix->solveSystem(matrix, b1, &opCount);

	delete[] weightMomentums;
	delete matrix;
	delete[] b1; 
	delete[] roots;

	return bigA;
	}
	delete[] weightMomentums;
	delete matrix;
	delete[] b1;
	throw "Problem with cubic roots!";
	return NULL;
}
double* dummyIntegral::buildThreePointNewton(double a, double b, double* nodes){
	double* weightMomentums = new double[3];
	weightMomentums[0] = mu[0](b) - mu[0](a);
	weightMomentums[1] = mu[1](b) - mu[1](a);
	weightMomentums[2] = mu[2](b) - mu[2](a);
	SillyMatrix* matrix = new SillyMatrix(3, 3);
	(*matrix)(0, 0) = 1.0;					(*matrix)(0, 1) = 1.0;							(*matrix)(0, 2) = 1.0;
	(*matrix)(1, 0) = nodes[0];				(*matrix)(1, 1) = nodes[1];						(*matrix)(1, 2) = nodes[2];
	(*matrix)(2, 0) = nodes[0]*nodes[0];	(*matrix)(2, 1) = nodes[1]*nodes[1];			(*matrix)(2, 2) = nodes[2]*nodes[2];
	int opCount = 0;
	double* bigA = matrix->solveSystem(matrix, weightMomentums, &opCount);

	delete[] weightMomentums;
	delete matrix;

	return bigA;
}
double dummyIntegral::integrateWithComposedQuadraticNewton(double a, double b, int numberOfSteps){
	double* newtonNodes = new double[3];
	double result = 0.0;
	double start = a;
	double* approximationMesh = this->buildApproximationMesh(a, b, numberOfSteps);
	double* As;
	for(int i = 0; i < numberOfSteps - 2; i+=2){
		newtonNodes[0] = approximationMesh[i];
		newtonNodes[1] = approximationMesh[i+1];
		newtonNodes[2] = approximationMesh[i+2];
		As = this->buildThreePointNewton(newtonNodes[0], newtonNodes[2], newtonNodes);
		if(As == NULL){
			delete[] approximationMesh;
			throw "Can't find big A's for Newton";
			return -1.0;
		}
		
		result += As[0]*fx(newtonNodes[0]);
		result += As[1]*fx(newtonNodes[1]);
		result += As[2]*fx(newtonNodes[2]);
		delete[] As;
		
	}
	delete[] approximationMesh;
	delete[] newtonNodes;
	return result;
	
}
double dummyIntegral::integrateWithComposedQuadraticGauss(double a, double b, int numberOfSteps){
	double result = 0.0;
	double* approximationMesh = this->buildApproximationMesh(a, b, numberOfSteps);
	double* As;
		for(int i = 0; i < numberOfSteps - 2; i+=2){
			try{
			As = this->buildThreePointGaussianQuadratic(approximationMesh[i], approximationMesh[i+2]);
			}catch(char* s){
				delete[] approximationMesh;
				throw s;
				return INT_MIN;
			}
			if(As == NULL){
				delete[] approximationMesh;
				throw "Can't find big A's for Gauss";
				return INT_MIN;
			}
			result += As[0]*fx(approximationMesh[i]);
			result += As[1]*fx(approximationMesh[i+1]);
			result += As[2]*fx(approximationMesh[i+2]);
			delete[] As;
		}
		delete[] approximationMesh;
	return result;
}

double dummyIntegral::integrateWithComposedQuadratic(double a, double b, Q_TYPE quadraticType){
	double result = 0.0;
		if(quadraticType == Q_NEWTON)
			result = this->integrateWithComposedQuadraticNewton(a, b, 61);
		else{
			try{
			result = this->integrateWithComposedQuadraticGauss(a, b, 111);
			}catch(char* s){
				cout <<"Exception: " << s << endl;
				return INT_MIN;
			}
		}
	return result;
}
typedef double (dummyIntegral::*integralCore)(double, double, int);
double dummyIntegral::integrateWithComposedQuadratic(double a, double b, double epsilon, Q_TYPE quadraticType){
	int numberOfStepsH1 = 11;
	int numberOfStepsH2 = 2*numberOfStepsH1 - 1;
	double m = (quadraticType == Q_NEWTON ? 7.0 : 63.0);
	double rungeDelta = 1;
	double S1, S2;
	integralCore integral = (quadraticType == Q_NEWTON ? &dummyIntegral::integrateWithComposedQuadraticNewton : &dummyIntegral::integrateWithComposedQuadraticGauss);
	while(rungeDelta > epsilon){
		try{
			S1 = (this->*integral)(a, b, numberOfStepsH1);
			S2 = (this->*integral)(a, b, numberOfStepsH2);
		}catch(char* s){
			cout <<"Exception: "<< s << endl;
			return INT_MIN;
		}
		rungeDelta = fabs((S2 - S1))/m;
		cout << "Delta: " << rungeDelta << endl;
		numberOfStepsH1 = numberOfStepsH2;
		numberOfStepsH2 *= 2;
		numberOfStepsH2 -= 1;
	}
	cout << numberOfStepsH2 << endl;
	return S2;
}
double dummyIntegral::processAitken(double a, double b, double epsilon, Q_TYPE quadraticType){
	integralCore integral = (quadraticType == Q_NEWTON ? &dummyIntegral::integrateWithComposedQuadraticNewton : &dummyIntegral::integrateWithComposedQuadraticGauss);
	int numberOfStepsH1 = 11;
	int numberOfStepsH2 = 2*numberOfStepsH1 - 1;
	int numberOfStepsH3 = 2*numberOfStepsH2 - 1;
	double delta = 0.0;
	double Cm = 1;
	double ape = 0.0;
	double S1, S2, S3;
	cout << (quadraticType == Q_NEWTON ? "Newton Aitken process simulation: " : "Gauss Aitken process simulation: ") << endl; 
	do{
	try{
		S1 = (this->*integral)(a, b, numberOfStepsH1);
		S2 = (this->*integral)(a, b, numberOfStepsH2);
		S3 = (this->*integral)(a, b, numberOfStepsH3);
	}catch(char* s){
		cout <<"Exception: " << s << endl;
		return (double)INT_MIN;
	}
	
	Cm = fabs((S1-S2)/(S2-S3));
	ape = log(Cm)/log(2.0);
	cout <<"APE: " <<  ape << endl;
	delta = Cm*pow((b-a)/numberOfStepsH3, ape);
	numberOfStepsH1 = numberOfStepsH3;
	numberOfStepsH2 = 2*numberOfStepsH1 -1;
	numberOfStepsH3 = 2*numberOfStepsH2 - 1;
	}while(delta > epsilon);
	cout << "end of Aitken process simulation" << endl;
	return ape;

}
double dummyIntegral::processAitken(double a, double b, Q_TYPE quadraticType){
	integralCore integral = (quadraticType == Q_NEWTON ? &dummyIntegral::integrateWithComposedQuadraticNewton : &dummyIntegral::integrateWithComposedQuadraticGauss);
	int numberOfStepsH1 = 11;
	int numberOfStepsH2 = 2*numberOfStepsH1 - 1;
	int numberOfStepsH3 = 2*numberOfStepsH2 - 1;
	double Cm = 1;
	double ape = 0.0;
	double theoreticalAPE = (quadraticType == Q_NEWTON ? 3.0 : 6.0);
	double S1, S2, S3;

	do{
	try{
		S1 = (this->*integral)(a, b, numberOfStepsH1);
		S2 = (this->*integral)(a, b, numberOfStepsH2);
		S3 = (this->*integral)(a, b, numberOfStepsH3);
	}catch(char* s){
		throw s;
		return (double)INT_MIN;
	}
	
	Cm = fabs((S1-S2)/(S2-S3));
	ape = log(Cm)/log(2.0);
	numberOfStepsH1 = numberOfStepsH3;
	numberOfStepsH2 = 2*numberOfStepsH1 -1;
	numberOfStepsH3 = 2*numberOfStepsH2 - 1;
	}while(ape < theoreticalAPE);
	return ape;
}
double dummyIntegral::integrateWithOptimalStepThroughAitken(double a, double b, double epsilon, Q_TYPE quadraticType){
	int numberOfStepsH1 = 11;
	int numberOfStepsH2 = 2*numberOfStepsH1 - 1;
	double Ho = 0;
	double m = 0;
	double rungeDelta = 1;
	double S1, S2;
	integralCore integral = (quadraticType == Q_NEWTON ? &dummyIntegral::integrateWithComposedQuadraticNewton : &dummyIntegral::integrateWithComposedQuadraticGauss);
		try{
			m = this->processAitken(a, b, quadraticType);
			S1 = (this->*integral)(a, b, numberOfStepsH1);
			S2 = (this->*integral)(a, b, numberOfStepsH2);
		}catch(char* s){
			cout <<"Exception: " << s << endl;
			return INT_MIN;
		}
		Ho = (b-a)/numberOfStepsH2*pow( ( epsilon*(1 - pow(2.0, -m) )/(fabs(S2 - S1))), 1.0/m);
	numberOfStepsH1 = (int)ceil(((b-a)/Ho));
	numberOfStepsH2 = 2*numberOfStepsH1 - 1;
	cout << "Optimal steps: " << numberOfStepsH1 << endl;
	cout << "Double the number of optimal steps: " << numberOfStepsH2 << endl;
	try{
			S1 = (this->*integral)(a, b, numberOfStepsH1);
			S2 = (this->*integral)(a, b, numberOfStepsH2);
		}catch(char* s){
			cout << s << endl;
			return INT_MIN;
		}
	cout << "With optimal number of steps: " << S1 << endl;
	cout << "With double the optimal steps: " << S2 << endl;
	return S2;
}
double dummyIntegral::integrateWithOptimalStep(double a, double b, double epsilon, Q_TYPE quadraticType){
	int numberOfStepsH1 = 11;
	int numberOfStepsH2 = 2*numberOfStepsH1 - 1;
	double Ho = 0;
	double m = 0;
	double rungeDelta = 1;
	double S1, S2;
	integralCore integral = (quadraticType == Q_NEWTON ? &dummyIntegral::integrateWithComposedQuadraticNewton : &dummyIntegral::integrateWithComposedQuadraticGauss);
		m = (quadraticType == Q_NEWTON ? 3.0 : 6.0);
		try{
			S1 = (this->*integral)(a, b, numberOfStepsH1);
			S2 = (this->*integral)(a, b, numberOfStepsH2);
		}catch(char* s){
			cout << s << endl;
			return INT_MIN;
		}
		Ho = (b-a)/numberOfStepsH2*pow( ( epsilon*(1 - pow(2.0, -m) )/(fabs(S2 - S1))), 1.0/m);
	numberOfStepsH1 = (int)ceil(((b-a)/Ho));

	numberOfStepsH2 = 2*numberOfStepsH1 - 1;
	cout << "Optimal steps: " << numberOfStepsH1 << endl;
	cout << "Double the number of optimal steps: " << numberOfStepsH2 << endl;
	try{
			S1 = (this->*integral)(a, b, numberOfStepsH1);
			S2 = (this->*integral)(a, b, numberOfStepsH2);
		}catch(char* s){
			cout <<"Exception: " << s << endl;
			return INT_MIN;
		}
	cout << "With optimal number of steps: " << S1 << endl;
	cout << "With double the optimal steps: " << S2 << endl;
	return S2;
}

int _tmain(int argc, _TCHAR* argv[])
{
	calc_func* mus = new calc_func[6];
	mus[0] = &i0_mu0;
	mus[1] = &i0_mu1;
	mus[2] = &i0_mu2;
	mus[3] = &i0_mu3;
	mus[4] = &i0_mu4;
	mus[5] = &i0_mu5;
	dummyIntegral* dummy = new dummyIntegral(&i0_fx, &i0_px, mus, 6);
	double res = 0.0;
	cout.precision(12);
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << endl;
	cout << "    3.3                                 " << endl;
	cout << "    /-\\ " << endl;
	cout << "    |" << endl;
	cout << "    |" << endl;
	cout << "    |" << endl;
	cout << "    |(x-2.1)^(-0.4)*(4.5cos(7x)*exp(-2x/3) + 1.4sin(1.5x)*exp(-x/3)+3)dx" << endl;
	cout << "    |" << endl;
	cout << "    |" << endl;
	cout << "    |" << endl;
	cout << "  \\-/ " << endl;
	cout << "   2.1                                      " << endl;
	cout << endl;
	cout << "12 digit precision from Wolfram: 4.46151270533" << endl;
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	try{
		res = dummy->integrateWithComposedQuadratic(2.1, 3.3, Q_NEWTON);
	cout <<"Integration with Newton Qudratic constant step result: "  << res << endl;
	}catch(char* s){
		cout <<"Exception: " << s << endl;
	}
	res = dummy->integrateWithComposedQuadratic(2.1, 3.3, Q_GAUSS);
	cout <<"Integration with Gauss Qudratic constant step result: " << res << endl;
	res = dummy->integrateWithComposedQuadratic(2.1, 3.3, 1e-12, Q_NEWTON);
	cout << "Integration with Runge rule usage (Newton): " << res << endl;
	res = dummy->integrateWithComposedQuadratic(2.1, 3.3, 1e-12, Q_GAUSS);
	cout << "Integration with Runge rule usage (Gauss): " << res << endl;
	dummy->processAitken(2.1, 3.3, 1e-12, Q_NEWTON);
	dummy->processAitken(2.1, 3.3, 1e-12, Q_GAUSS);
	dummy->integrateWithOptimalStep(2.1, 3.3, 1e-12, Q_NEWTON);
	//dummy->integrateWithOptimalStep(2.1, 3.3, 1e-12, Q_GAUSS);
	dummy->integrateWithOptimalStepThroughAitken(2.1, 3.3, 1e-12, Q_NEWTON);
	//dummy->integrateWithOptimalStepThroughAitken(2.1, 3.3, 1e-12, Q_GAUSS);
	delete mus;
	delete dummy;
	cin.get();
	return 0;
}

