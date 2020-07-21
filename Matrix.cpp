#include<iostream>
#include<math.h>
#include"Matrix.h"

using namespace std;

//class Matrix;
//double arrMult(Matrix A, Matrix B,  int row, int col);
//double det(Matrix A);





	Matrix::Matrix():row(0), col(0){}

	Matrix::Matrix(int rowSize, int colSize, double num):row(rowSize), col(colSize){
		arr =  new double*[row];
		for(int i=0; i<row; i++){
			arr[i] = new double[col];
		}
		for(int i=0; i<row; i++){
			for(int j=0; j<col; j++){
				arr[i][j] = num;
			}
		}
	}

	Matrix::Matrix(double **vals, int rowSize, int colSize):row(rowSize), col(colSize){
		arr = new double*[row];
		for(int i=0; i<row; i++){
			arr[i] = new double[col];
			for(int j=0; j<col; j++){
				arr[i][j] = vals[i][j];
			}
		}
	}

	Matrix::Matrix(const Matrix& copy):row(copy.row), col(copy.col){
		arr = new double*[row];
		for(int i=0; i<row; i++){
			arr[i] = new double[col];
		}
		for(int i=0; i<row; i++){
			for(int j=0; j<col; j++){
				arr[i][j] = copy.arr[i][j];
			}
		}
	}

	Matrix Matrix::T(){
		Matrix temp(col, row);
		for(int i=0; i<temp.row; i++){
			for(int j=0; j<temp.col; j++){
				temp.arr[i][j] = arr[j][i];
	//				cout<<i<<' '<<j<<endl;
			}
		}
		return temp;
	}


	Matrix::Matrix(const Matrix &A, const Matrix &B):row(A.row), col(B.col){
		arr = new double*[row];
		for(int i=0; i<row; i++){
			arr[i] = new double[col];
		}
		for(int i=0; i<row; i++){
//			cout<<"Row: "<<i<<endl;
			for(int j=0; j<col; j++){
				arr[i][j] = arrMult(A, B, i, j);
			}
		}
	}

	Matrix::Matrix(double num, Matrix &A):row(A.row), col(A.col){
		arr = new double*[row];
		for(int i=0; i<row; i++){
			arr[i] = new double[col];
		}
		for(int i=0; i<row; i++){
			for(int j=0; j<col; j++){
				arr[i][j] = num*A.arr[i][j];
			}
		}
	}

	Matrix::Matrix(Matrix A, int row_d, int col_d){
		row = A.row-1;
		col = A.col-1;
		arr = new double*[row];
		for(int i=0; i<row; i++){
			arr[i] = new double[col];
		}
//		int x=0, y=0;
		for(int i=0, x=0; i<row; i++, x++){
			if(x==row_d){
				x++;
			}
			for(int j=0, y=0; j<col; j++, y++){
				arr[i][j] = A.arr[x][y==col_d ? ++y : y];
			}
		}
	}


	Matrix::~Matrix(){
		for(int i=0; i<row; i++){
			delete[] arr[i];
		}
		delete[] arr;
	}


	void Matrix::operator=(const Matrix &copy){
		if(row!=0 && col!=0){
			for(int i=0; i<row; i++){
				delete[] arr[i];
			}
			delete[] arr;
		}
		row = copy.row;
		col = copy.col;
		arr = new double*[row];
		for(int i=0; i<row; i++){
			arr[i] = new double[col];
		}
		for(int i=0; i<row; i++){
			for(int j=0; j<col; j++){
				arr[i][j] = copy.arr[i][j];
			}
		}
	}

ostream& operator<<(ostream& os, const Matrix &matrix){
//	os << '['<<endl;
	for(int i=0; i<matrix.row; i++){
		for(int j=0; j<matrix.col; j++){
			os << matrix.arr[i][j]<<'\t';
		}
		os << endl;
	}
//	os<<']'<<endl;
	return os;
}

double arrMult(const Matrix &A, const Matrix &B, const int row, const int col){
	double temp=0;
	for(int i=0; i<A.col; i++){
		temp += A.arr[row][i] * B.arr[i][col];
	}
	return temp;
}

Matrix operator *(double const num, Matrix &A){
	Matrix temp(A);
	for(int i=0; i<temp.row; i++){
		for(int j=0; j<temp.col; j++){
			temp.arr[i][j] = num * A.arr[i][j];
		}
	}
	return temp;
}

Matrix operator *(Matrix A, Matrix B){
	Matrix temp(A,B);
	return temp;
}

Matrix operator+(Matrix A, Matrix B){
	Matrix temp(A);
	for(int i=0; i<temp.row; i++){
		for(int j=0; j<temp.col; j++){
			temp.arr[i][j] += B.arr[i][j];
		}
	}
	return temp;
}

Matrix operator-(Matrix A, Matrix B){
	Matrix temp(A);
	for(int i=0; i<temp.row; i++){
		for(int j=0; j<temp.col; j++){
			temp.arr[i][j] -= B.arr[i][j];
		}
	}
	return temp;
}

double det(Matrix A){
		if(A.row==2 && A.col==2){//only works on square matrices ofc
			return A.arr[0][0]*A.arr[1][1] - A.arr[1][0]*A.arr[0][1];
		}
		else{
			double ret=0;
			for(int i=0; i<A.col; i++){
				Matrix temp(A, 0, i);
				if(i%2==0){
					ret += A.arr[0][i] * det(temp);
				}
				else{
					ret -= A.arr[0][i] * det(temp);
				}
			}
			return ret;
		}
}

void setValue(double* A, double* B, int size){
	for(int i=0; i<size; i++){
		A[i] = B[i];
	}
}

int pow(int A, int B){
	if(B==0) return 1;
	else return A*pow(A, B-1);
}

Matrix inv(Matrix A){
	Matrix temp (A.col, A.row);
	double inv_det = det(A);
//	cout<< inv_det<<endl;
	if(A.row==2 && A.col==2){
		temp.arr[0][0] = 1.0/inv_det * A.arr[1][1];
		temp.arr[0][1] = -1.0/inv_det * A.arr[1][0];
		temp.arr[1][0] = -1.0/inv_det * A.arr[0][1];
		temp.arr[1][1] = 1.0/inv_det * A.arr[0][0];
	}
	else{
		Matrix C(A);
		for(int i=0; i<C.row; i++){
			for(int j=0; j<C.col; j++){
				Matrix sub(A,i,j);
				C.arr[i][j] = pow(-1,i+j) * det(sub);
				temp.arr[j][i] = C.arr[i][j] /inv_det;
			}
		}
		
	}
	return temp;
}




double case3(Matrix x, Matrix mu, Matrix Sigma, double Pw){
	Matrix Wi(1,1);
	Matrix wi(1,1);
	Matrix wi0(1,1);
	double wi0_cont;

	Matrix inv_Sigma = inv(Sigma);
	double det_Sigma = det(Sigma);

	Wi = ((x.T() * inv_Sigma)*x);
	Wi = -0.5 * Wi;

	wi = (inv_Sigma * mu).T() * x;

	wi0 = ((mu.T() * inv_Sigma)*mu);
	wi0 = -0.5 * wi0;

	wi0_cont = -0.5 * log(det_Sigma) + log(Pw);

	return Wi.arr[0][0] + wi.arr[0][0] + wi0.arr[0][0] + wi0_cont;

//	Wi = ((x.T() * inv(Sigma))*x);
//	Wi = -0.5 * Wi;

//	wi = (inv(Sigma) * mu).T() * x;

//	wi0 = ((mu.T() * inv(Sigma))*mu);
//	wi0 = -0.5 * wi0;

//	wi0_cont = -0.5 * log(det(Sigma)) + log(Pw);

//	return Wi.arr[0][0] + wi.arr[0][0] + wi0.arr[0][0] + wi0_cont;

//	double val = (-0.5*((x.T() * inv(Sigma))*x) + ((Sigma*mu).T() * x) - 0.5*((mu.T() * inv(Sigma))*mu)).arr[0][0]
//	return  val -0.5*log(det(Sigma)) + log(Pw);
}

void Identity(Matrix &A, double scalar){
	A = 0 * A;
	for(int i=0; i<A.row; i++){
		A.arr[i][i] = scalar;
	}
}






