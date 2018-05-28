#include<iostream>
#include<cmath>
#include<stdio.h>

#define N 100
#define a 0.5
#define Epsilon 1

using namespace std;
double ** A;
double *B;
double h = 1.0 / N;
double *Y;

void Jacobi()
{
	double *X_now = new double[N-1];
	double *X_next = new double[N-1];
	for ( int i = 0; i < N-1; i++ )
		X_now[i] = 0;
	int num = 0;
	while (num < 5000)
	{
		for ( int i = 0; i < N-1; i++ )
		{
			double sum = 0.0;
			for ( int j = 0; j < N-1; j++ )
				if ( j != i )
					sum += A[i][j] * X_now[j]; 
			X_next[i] = (B[i] - sum)/A[i][i];
		}
		for ( int i = 0; i < N-1; i++ )
			X_now[i] = X_next[i];
		
		num++;		
	}
	double sum = 0.0;
	for ( int i = 0; i < N-1; i++)
		sum += (X_now[i] - Y[i])*(X_now[i] - Y[i]);	
	sum = sqrt(sum);
	cout << "Jacobi迭代：" << endl;
	cout <<	"误差：" << sum << endl;
	for ( int i = 0; i < N-1; i++)
		printf("%.4f ",X_now[i]);
	cout << endl;
}

void Gauss_Seidel()
{
	double *X_now = new double[N-1];
	double *X_next = new double[N-1];
	for ( int i = 0; i < N-1; i++ )
		X_now[i] = 0;
	int num = 0;
	while (num < 5000)
	{
		for ( int i = 0; i < N-1; i++ )
		{
			double sum = 0.0;
			for ( int j = 0; j < i; j++ )
				sum += A[i][j] * X_next[j];
			for ( int j = i; j < N-1; j++ )
				sum += A[i][j] * X_now[j]; 
			X_next[i] = X_now[i] + ((B[i] - sum)/A[i][i]);
		}
		for ( int i = 0; i < N-1; i++ )
			X_now[i] = X_next[i];
		num++;		
	}
	double sum = 0.0;
	for ( int i = 0; i < N-1; i++)
		sum += (X_now[i] - Y[i])*(X_now[i] - Y[i]);	
	sum = sqrt(sum);
	cout << "Gauss_Seidel迭代：" << endl;
	cout <<	"误差：" << sum << endl;
	for ( int i = 0; i < N-1; i++)
		printf("%.4f ",X_now[i]);
	cout << endl;
}

void SOR()
{
		double *X_now = new double[N-1];
	double *X_next = new double[N-1];
	for ( int i = 0; i < N-1; i++ )
		X_now[i] = 0;
	int num = 0;
	while (num < 5000)
	{
		for ( int i = 0; i < N-1; i++ )
		{
			double sum = 0.0;
			for ( int j = 0; j < i; j++ )
				sum += A[i][j] * X_next[j];
			for ( int j = i; j < N-1; j++ )
				sum += A[i][j] * X_now[j]; 
			X_next[i] = X_now[i] + 1.5*((B[i] - sum)/A[i][i]);
		}
		for ( int i = 0; i < N-1; i++ )
			X_now[i] = X_next[i];
		num++;		
	}
	double sum = 0.0;
	for ( int i = 0; i < N-1; i++)
		sum += (X_now[i] - Y[i])*(X_now[i] - Y[i]);	
	sum = sqrt(sum);
	cout << "SOR迭代：" << endl;	
	cout << "误差：" << sum << endl;
	for ( int i = 0; i < N-1; i++)
		printf("%.4f ",X_now[i]);
	cout << endl;
}

int main()
{
	Y = new double[N-1];
	for ( int i = 0; i < N-1; i++ )
		Y[i] = (1-a)/(1-exp(-1/Epsilon)) *(1-exp(-((i+1) * h)/Epsilon)) + a*((i+1) * h);
	cout << "精确解："<< endl; 
	for ( int i = 0; i < N-1; i++)
		printf("%.4f ",Y[i]);
	cout << endl;	
	A = new double* [N-1];
	for ( int i = 0; i < N-1; i++ )
		A[i] = new double[N-1];
	for ( int i = 0; i < N-1; i++)
		for (int j = 0; j < N-1; j++ )
			A[i][j] = 0.0;
			
	for (int i = 0; i < N-1; i++ )
		A[i][i] = -(2.0*Epsilon+h);
	for (int i = 1; i < N-1; i++ )
		A[i][i-1] = Epsilon;
	for (int i = 0; i < N-2; i++ )
		A[i][i+1] = Epsilon+h;
//	for ( int i = 0; i < N-1; i++)
//		{
//			for (int j = 0; j < N-1; j++ )
//				cout << A[i][j] << ' ';
//			cout << endl;
//		}
	B = new double[N-1];
	for ( int i = 0; i < N-1; i++ )
		B[i] = a * h * h;
	B[N-2] -= (Epsilon+h); 
//	for ( int i = 0; i < N-1; i++ )
//		cout << B[i] << ' ';
	Jacobi();
	Gauss_Seidel();
	SOR();
		
	return 0;
}
