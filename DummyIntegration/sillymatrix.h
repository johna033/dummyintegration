#include "stdafx.h"
#include <vector>
#include <iostream>
using namespace std;

#define L_POS 0
#define U_POS 1
#define RP_POS 2
#define ORIG_POS 3
#define EPSILON 1e-10

class SillyMatrix{
private:
	int cols, rows;
	double** matrix;

public: 
	SillyMatrix(int cols, int rows):cols(cols), rows(rows){
		matrix = new double*[rows];
		for(int i = 0; i < rows; i++)
			matrix[i] = new double[cols];

	};
	SillyMatrix(SillyMatrix* A): cols(A->cols), rows(A->rows){
		matrix = new double*[rows];
		for(int i = 0; i < rows; i++){
			matrix[i] = new double[cols];
			for(int j = 0; j < cols; j++)
				matrix[i][j] = (*A)(i, j);
		}
	}
	SillyMatrix():cols(0), rows(0){};
	~SillyMatrix();
	SillyMatrix* operator*(SillyMatrix* a);
	double& operator()(int row, int col){return matrix[row][col];}
	vector<SillyMatrix*>* LUPDecomposition(int* toggles, int* opCount);
	SillyMatrix* getInverse();
	double* solveSystem(SillyMatrix* A, double* b, int *opCount);
	double* solveSystemWithLU(SillyMatrix* LU, SillyMatrix* P, double* b);
	void printOut();
	double frobeniusNorm();
	double determinant();
	int getDim();
	int getRank();
private:
	SillyMatrix(SillyMatrix& A){};
	void swapColumns(int i, int j);
	void swapRows(int i, int j);
	SillyMatrix* identityMatrix();
	SillyMatrix* zeroMatrix();
	vector<SillyMatrix*>* sliceUpUPResult(SillyMatrix* result);
	double* solve(SillyMatrix* LU, double* b);
	double* multiplyMatrixByVector(SillyMatrix* A, double* vec);
	bool isRowZero(int j, SillyMatrix* A);
	void deleteVector(vector<SillyMatrix*>* vec);
};