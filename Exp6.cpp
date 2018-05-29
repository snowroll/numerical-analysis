#include <iostream>
#include <cmath>
#include <stdio.h>
using namespace std;

#define Epsilon 1
#define a 0.5
#define N 100

double h = 1.0 / N;
double **A;  //系数矩阵
double *B;  //结果矩阵
double *Y;  //精确解

void Deviation(double *X, double *Y, char method[]){
    double dev = 0.0;
    cout << method << " result: " << endl;
    for(int i = 0; i < N - 1; i++){
        printf("%.4f  *",X[i]);  //输出解
        dev += (X[i] - Y[i]) * (X[i] - Y[i]);
    }
    dev = sqrt(dev);
    cout << endl;
    cout << method << " deviation: " << dev << endl;
}

void Jacobi(double *Y){
    double *X_k = new double[N - 1];
    double *X_k_1 = new double[N - 1];
    for(int i = 0; i < N - 1; i++)  //初始化迭代的X
        X_k[i] = 0.0;
    int iter_num = 0;
    while(iter_num <= 10000){
        for(int i = 0; i < N - 1; i++){
            double tmp = 0.0;
            for(int j = 0; j < N - 1; j++)
                if(i != j)
                    tmp -= A[i][j] * X_k[j];
            X_k_1[i] = (tmp + B[i]) / A[i][i];
        }
        for(int i = 0; i < N -1; i++)
            X_k[i] = X_k_1[i];
        iter_num++;
    }
    char method[] = "Jacobi";
    Deviation(X_k, Y, method);
}

void Gauss_Seidel(double *Y){
    double *X_k = new double[N - 1];
    double *X_k_1 = new double[N - 1];
    for(int i = 0; i < N - 1; i++)  //初始化迭代的X
        X_k[i] = 0.0;
    int iter_num = 0;
    while(iter_num <= 10000){
        for(int i = 0; i < N - 1; i++){
            double tmp = 0.0;
            for(int j = 0; j < i; j++)
                tmp -= A[i][j] * X_k_1[j];
            for(int j = i + 1; j < N - 1; j++)
                tmp -= A[i][j] * X_k[j];
            X_k_1[i] = (B[i] + tmp) / A[i][i];  
        }
        for(int i = 0; i < N -1; i++)
            X_k[i] = X_k_1[i];
        iter_num++;
    }
    char method[] = "Gauss_Seidel";
    Deviation(X_k, Y, method);
}

void  SOR(double *Y){
    double *X_k = new double[N - 1];
    double *X_k_1 = new double[N - 1];
    for(int i = 0; i < N - 1; i++)  //初始化迭代的X
        X_k[i] = 0.0;
    int iter_num = 0;
    while(iter_num <= 10000){
        for(int i = 0; i < N - 1; i++){
            double tmp = 0.0;
            for(int j = 0; j < i; j++)
                tmp -= A[i][j] * X_k_1[j];
            for(int j = i; j < N - 1; j++)
                tmp -= A[i][j] * X_k[j];
            X_k_1[i] = X_k[i] + (B[i] + tmp) / A[i][i];  
        }
        for(int i = 0; i < N -1; i++)
            X_k[i] = X_k_1[i];
        iter_num++;
    }
    char method[] = "SOR";
    Deviation(X_k, Y, method);
}

int main(){
    A = new double*[N - 1];
    double const_num[3] = {Epsilon + h, -h - 2.0 * Epsilon, Epsilon};
    for(int i = 0; i < N - 1; i++){  //初始化A矩阵
        A[i] = new double[N - 1];
        for(int j = 0; j < N -1; j++)
            if((i - j) >= -1 && (i - j) <= 1)
                A[i][j] = const_num[i-j+1];
            else 
                A[i][j] = 0.0;
    }

    B = new double[N - 1];
    for(int i = 0; i < N - 1; i++)
        B[i] = a * h * h;
    B[N - 2] -= (h + Epsilon);  

    Y = new double[N - 1];
    for(int i = 0; i < N - 1; i++)
        Y[i] = (1 - a) / (1 - exp(-1 / Epsilon)) * (1 - exp(-((i + 1) * h) / Epsilon)) + a * (i + 1) * h;
    Jacobi(Y);
    Gauss_Seidel(Y);
    SOR(Y);
    return 0;
}