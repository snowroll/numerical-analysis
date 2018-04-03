#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

#define N 7
#define M 4
const int num = 200;  //样例点个数

float x[N] = {-1, -0.5, 0, 0.5, 1, 1.5, 2}, 
      y[N] = {-4.467, -0.452, 0.551, 0.048, -0.447, 0.549, 4.552},
      phi[M][M] = {0},
      d[M] = {0},
      coef[M] = {0};  //解的存储数组

void gaussin(float a[M][M], float b[]){  //高斯消元求解方程组
    for(int i = 0; i < M; i++)  //对角元素不能为0
        if (a[i][i] == 0){
            cout << "can't use gaussin meathod" << endl;
            return;
        }

    float m[M];  //存储初等行变换的系数，用于行的相减
    for(int k = 0; k < M - 1; k++){  //化成对角矩阵
        for(int i = k + 1; i < M; i++)  
            m[i] = a[i][k] / a[k][k];

        for(int i = k + 1; i < M; i++){ 
            for(int j = 0; j < M; j++){
                a[i][j] = a[i][j] - m[i] * a[k][j];
            }
            b[i] = b[i] - m[i] * b[k];
        }
    }

    coef[M - 1] = b[M - 1] / a[M - 1][M - 1];  //先计算出最后一个未知数
    for(int i = M - 2; i >= 0; i--){  //求出每个未知数的值
        float sum = 0;
        for(int j = i + 1; j < M; j++)
            sum += a[i][j] * coef[j];
        coef[i] = (b[i] - sum) / a[i][i];
    }
}

float Expre(float _x){
	float ans = 0;
	switch(M){
		case 2:
			ans = coef[0] + coef[1] * _x;
			break;
		case 3:
			ans = coef[0] + coef[1] * _x + coef[2] * _x * _x;
			break;
		case 4:
			ans = coef[0] + coef[1] * _x + coef[2] * _x * _x + coef[3] * _x * _x * _x;
			break;
	}
	return ans;
}

int main(){
	for(int i = 0; i < M; i++)
		for(int j = 0; j < M; j++)
			for(int k = 0; k < N; k++)
				phi[i][j] += pow(x[k], i) * pow(x[k], j);
			
	for(int i = 0; i < M; i++)
		for(int k = 0; k < N; k++)
			d[i] += pow(x[k], i) * y[k];

	gaussin(phi, d);
	
	if(M == 2){ 
		cout << "y = " << coef[0] << " + " << coef[1] << "x" << endl;
		cout << "error: ";
		for(int i = 0; i < N; i++)
			cout << Expre(x[i]) - y[i] << " ";
		cout << endl;
	}
	else if(M == 3){
		cout << "y = " << coef[0] << " + " << coef[1] << "x" << " + " << coef[2] << "x2" << endl;
		cout << "error: ";
		for(int i = 0; i < N; i++)
			cout << Expre(x[i]) - y[i] << " ";
		cout << endl;
	}
	else if(M == 4){
		cout << "y = " << coef[0] << " + " << coef[1] << "x" << " + " << coef[2] << "x2" << " + " << coef[3] << "x3" << endl;
		cout << "error: ";
		for(int i = 0; i < N; i++)
			cout << Expre(x[i]) - y[i] << " ";
		cout << endl;
	}


	float h = (x[N-1] - x[0]) / num;  //step length
	float _x[num + 1], _y[num + 1];
	ofstream outx("square-x.txt");
	ofstream outy("square-y.txt");
	if(outx.is_open() && outy.is_open()){ 
		for(int i = 0; i < num; i++){
			_x[i] = x[0] + i * h;
			_y[i] = Expre(_x[i]);
			outx << _x[i] << endl;
			outy << _y[i] << endl;
		}
	}
	outx.close();
	outy.close();
	
	return 0;
}