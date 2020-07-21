//Matrix h

#ifndef MATRIX_H
#define MATRIX_H

#include<iostream>

class Matrix{
	public:
		int row, col;
		double **arr;

		Matrix();//generates an empty 0x0 matrix
		Matrix(int rowSize, int colSize, double num=0);//generates a row x col matrix, all values = num.  default num=0
		Matrix(const Matrix& copy);//generates a copy-Matrix
		Matrix(double **vals, int rowSize, int colSize);

		Matrix T();//generates transpose of matrix
		Matrix(const Matrix &A, const Matrix &B);//generates matrix of A*B
		Matrix(double num, Matrix &A);//generates scalar num * Matrix A

		Matrix(Matrix A, int row_d, int col_d);//generates Matrix A, excluding given row and col.  used for determinant and inverse
		
		~Matrix();

		void operator=(const Matrix &copy);//sets Matrix = copy
};

std::ostream& operator<<(std::ostream &os, const Matrix &matrix);//printing Matrix

Matrix operator == (Matrix A, Matrix B);

Matrix operator * (double const num, Matrix &A);//scalar num * Matrix A
Matrix operator * (Matrix A, Matrix B);//Matrix A * Matrix B

Matrix operator + (Matrix A, Matrix B);//Matrix A + Matrix B.  Must be legal computation.
Matrix operator - (Matrix A, Matrix B);//Matrix A - Matrix B.  Must be legal computation

void Identity(Matrix &A, double scalar);

double det(Matrix A);//returns determinant of legal Matrix.

Matrix inv(Matrix A);//returns inverse of legal Matrix.

double arrMult(const Matrix &A, const Matrix &B, const int row, const int col);//helper function for finding det and inv

class Vector{
	public:
		Matrix M;
		Vector();//generates a 2x1 Matrix.  Used for holding a coordinate vector
};

void setValue(double *A, double *B, int size);//sets 1D array A to 1D array B
int pow(int A, int B);//returns A^B

double case3(Matrix x, Matrix mu, Matrix Sigma, double Pw);//returns Case 3 discriminant value given x, mu, Sigma, and Pw


#endif
