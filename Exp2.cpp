#include <iostream>
#include <fstream>
using namespace std;

#define N 20  //插值的次数

float wn_1(int n, float array[], float x){
	float ans = 1.0; 
	for(int i = 0; i <= n; i++)
		ans *= (x - array[i]);
	return ans;
}

float _wn_1(int k, int num, float array[]){
	float ans = 1;
	for(int i = 0; i <= num; i++){
		if(i != k)
			ans *= (array[k] - array[i]);
	}
	return ans;
}

float Ln_x(float x[], float y[], int num, float _x){  //lagrange插值
	float Ln_x = 0;
	float wn_1_ans = wn_1(num, x, _x);
	for(int i = 0; i <= num; i++){
		if(_x - x[i] == 0){
			Ln_x = y[i];
			break;
		}
		Ln_x += y[i] * wn_1_ans / ((_x - x[i]) * _wn_1(i, num, x));
	}
	return Ln_x;
}

int main(){
	float h = 10.0 / N;  //step length
	float x[N + 1];
	float y[N + 1];
	float _x;
	for(int i = 0; i <= N; i++){
		x[i] = - 5.0 + i * h;
		y[i] = 1.0 / (1.0 + 16 * x[i] * x[i]);
	}
	ofstream outx("Lagrange-x.txt");
	ofstream outy("Lagrange-y.txt");
	if(outx.is_open() && outy.is_open()){
		int m = 1000;
		float len = 10.0 / m; float cur_x = -5.0 - len;
		for(int i = 0; i <= m; i++){
			cur_x += len; 
			outx << cur_x << endl;
			outy << Ln_x(x, y, N, cur_x) << endl;
		}
	}
	outx.close();
	outy.close();

	return 0;
}