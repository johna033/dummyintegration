#include "stdafx.h"
#include "sillymatrix.h"

SillyMatrix::~SillyMatrix(){
	if(rows*cols != 0){
		for(int i = 0; i < rows; i++)
			delete[] matrix[i];
		delete[] matrix;
	}
}
double SillyMatrix::frobeniusNorm(){
	double norm = 0.0;
	for(int i = 0; i < this->rows; i++){
		for(int j = 0; j < this->cols; j++){
			norm+=this->matrix[i][j]*this->matrix[i][j];
			}
		}
	return sqrt(norm);
}
int SillyMatrix::getDim(){
	return rows;
}
bool SillyMatrix::isRowZero(int j, SillyMatrix* A){
	int n = A->getDim();
	for(int i = 0; i < n; i++)
		if((*A)(j, i) != 0)
			return false;
	return true;
}
double* SillyMatrix::solveSystem(SillyMatrix* A, double* b, int* opCount){
	int toggles;
	vector<SillyMatrix*>* lup = A->LUPDecomposition(&toggles, opCount);
	SillyMatrix* P = (*lup)[RP_POS];
	SillyMatrix* LU = (*lup)[ORIG_POS];
	double* pb = this->multiplyMatrixByVector(P, b);
	double* x = NULL;
	 x = this->solve(LU, pb);
	 this->deleteVector(lup);
	return x;
}
double* SillyMatrix::solveSystemWithLU(SillyMatrix* LU, SillyMatrix* P, double* b){
	double* pb = this->multiplyMatrixByVector(P, b);
	double* x = NULL;
	 x = this->solve(LU, pb);
	 delete[] pb;
	return x;
}
vector<SillyMatrix*>* SillyMatrix::sliceUpUPResult(SillyMatrix* result){
	SillyMatrix* L = this->identityMatrix();
	SillyMatrix* U = this->zeroMatrix();
	for(int i = 0; i < this->rows; i++){
		for(int j = 0; j < this->cols; j++){
		if(i > j){
			(*L)(i, j) = (*result)(i, j);
		}
		if(i <= j){
			(*U)(i, j) = (*result)(i, j);
		}
		}
	}
	vector<SillyMatrix*>* cRes = new vector<SillyMatrix*>();
	cRes->push_back(L);
	cRes->push_back(U);
	return cRes;

}
double* SillyMatrix::multiplyMatrixByVector(SillyMatrix* A, double* vec){
	int n = A->getDim();
	double* result = new double[n];
	for(int i = 0; i < n; i++){
		result[i] = 0;
		for(int j = 0; j < n; j++){
			result[i] += vec[j]*(*A)(i, j);
		}
	}
	return result;
}
double SillyMatrix::determinant(){
	int rank, toggles;
	double res = 1.000;
	vector<SillyMatrix*>* lup = this->LUPDecomposition(&rank, &toggles);
	SillyMatrix* U = (*lup)[U_POS];
	for(int i = 0; i < this->rows; i++){
			res *= (*U)(i, i);
	}
	if(toggles % 2 != 0)
		res *= -1.0;
	this->deleteVector(lup);
	return res;

}
void SillyMatrix::deleteVector(vector<SillyMatrix*>* vec){
	int sz = vec->size();
	for(int i = 0; i < sz; i++){
		delete (*vec)[i];
		(*vec)[i] = NULL;
	}
	delete vec;
	vec = NULL;
}
double* SillyMatrix::solve(SillyMatrix* LU, double* b){
	int n = LU->getDim();
	double* x = new double[n];
	bool isBZero = true;
	for(int i = 0; i < n; i++){
		x[i] = b[i];
		if(x[i] != 0)
			isBZero = false;
	}
	if(isBZero)
		return x;
	  for (int i = 1; i < n; ++i)
		{
			double sum = x[i];
				for (int j = 0; j < i; ++j)
					sum -= (*LU)(i, j) * x[j];
			x[i] = sum;
		}
		x[n - 1] /= (*LU)(n-1, n-1);
		for (int i = n - 2; i >= 0; --i)
		 {
			double sum = x[i];
			for (int j = i + 1; j < n; ++j)
			sum -= (*LU)(i, j) * x[j];
			if(fabs((*LU)(i, i)) <= EPSILON){
			return NULL;
			}
			x[i] = sum / (*LU)(i, i);
			
		}
		return x;
	}
SillyMatrix* SillyMatrix::getInverse(){
	int rank, toggles;
	vector<SillyMatrix*>* lup = this->LUPDecomposition(&rank, &toggles);
	SillyMatrix* P = (*lup)[RP_POS];
	SillyMatrix* LU = (*lup)[ORIG_POS];
	SillyMatrix* result = new SillyMatrix(this->cols, this->rows);
	double* tmp = new double[this->rows];
	for(int i = 0; i < this-> cols; i++){
		for(int j = 0; j < this->rows; j++){
			tmp[j] = (*P)(j, i);
		}
		double* sol = this->solve(LU, tmp);
		if(sol == NULL)
			return result;
			for(int w = 0; w < this->rows; w++){
				(*result)(w, i) = sol[w];
		}
		delete[] sol;
	}
	delete[] tmp;
	return result;
}
int SillyMatrix::getRank(){
	SillyMatrix rankRevealing = (*this);
	int virtualPivotRow = 0;
	int virtualPivotCol = 0;
	int n = this->rows;
	for(int i = 0; i < n; i++){

			virtualPivotRow = i;
			virtualPivotCol = i;
			for(int col = i; col < n; col++){
			for( int row = i; row < n; row++ ) {
            if( fabs(rankRevealing( row,  col )) > EPSILON ) {
                virtualPivotRow = row;
				virtualPivotCol = col;
				break;
            }
			}
			}
			rankRevealing.swapColumns(virtualPivotCol, i);
			rankRevealing.swapRows(virtualPivotRow, i);
			for(int j = i+1; j < n; j++) {
            rankRevealing( j, i ) /= rankRevealing(i, i );
            for(int k = i+1; k < n; k++) {
                rankRevealing(j, k ) -= (rankRevealing( j, i) * rankRevealing(i, k));
			}
        }
	}
	int rank = 0;
	vector<SillyMatrix*>* rankRes = this->sliceUpUPResult(&rankRevealing);
	for(int i = 0; i < n; i++)
		if(!this->isRowZero(i, (*rankRes)[U_POS]))
			rank++;
	return rank;

}
vector<SillyMatrix*>* SillyMatrix::LUPDecomposition(int* toggles, int* opCount){
	SillyMatrix* result = new SillyMatrix(this);
	const int n = this->rows;
	SillyMatrix* e = this->identityMatrix();
	vector<SillyMatrix*>* completeResult;
	(*toggles) = 0;
    for( int i = 0; i < n; i++ ) {
        double pivotValue = 0;
        int pivotRow = i;
		
		for(int row = i; row < n; row++){
            if( fabs((*result)( row,  i )) > pivotValue ) {
                pivotValue = fabs((*result)(row,  i ));
                pivotRow = row;
        }
		}
		if(pivotRow != i){
			e->swapRows(pivotRow, i);
			result->swapRows(pivotRow, i);
			(*toggles)++;
		}
        if(pivotValue > EPSILON){
        for(int j = i+1; j < n; j++) {
            (*result)( j, i ) /= (*result)(i, i );
			(*opCount)++;
            for(int k = i+1; k < n; k++) {
                (*result)(j, k ) -= ((*result)( j, i) * (*result)(i, k));
				(*opCount) += 2; 
			}
        }
		}
    }
	completeResult = this->sliceUpUPResult(result);
	completeResult->push_back(e);
	completeResult->push_back(result);
	return completeResult;
}
SillyMatrix* SillyMatrix::identityMatrix(){
	SillyMatrix* result = new SillyMatrix(this->cols, this->rows);
	for(int i = 0; i < this->rows; i++){
		for(int j = 0; j < this->cols; j++){
			if(i == j){
				(*result)(i, j) = 1;
			}else{
				(*result)(i, j) = 0;
			}
		}
	}
	return result;
}
SillyMatrix* SillyMatrix::zeroMatrix(){
	SillyMatrix* result = new SillyMatrix(this->cols, this->rows);
	for(int i = 0; i < this->rows; i++){
		for(int j = 0; j < this->cols; j++){
			(*result)(i, j) = 0;
		}
	}
	return result;
}
SillyMatrix* SillyMatrix::operator*(SillyMatrix* a){
	SillyMatrix* result = new SillyMatrix(this->cols, this->rows);
	for(int i = 0;i < this->rows;i++){
		 for(int j = 0;j < this->cols;j++){
               (*result)(i,j) = 0.0;
			   for(int k=0;k < this->rows;k++){
               (*result)(i, j) = (*result)(i, j) + this->matrix[i][k] * (*a)(k, j);
          }
          }
     }
	return result;
}

void SillyMatrix::swapColumns(int i, int j){
	for(int x = 0; x < this->rows; x++){
		double t = this->matrix[x][i];
		this->matrix[x][i] = this->matrix[x][j];
		this->matrix[x][j] = t;
	}
}
void SillyMatrix::swapRows(int i, int j){
	for(int x = 0; x < this->cols; x++){
		double t = this->matrix[i][x];
		this->matrix[i][x] = this->matrix[j][x];
		this->matrix[j][x] = t;
	}
}
void SillyMatrix::printOut(){
	char result[100]; 
	for(int i = 0; i < this->rows; i++)
	{
		for(int j = 0; j < this->cols; j++){
			sprintf(result, "%2.3f", this->matrix[i][j]);
			cout << result << "  ";		
		}
		cout << endl;
	}
}