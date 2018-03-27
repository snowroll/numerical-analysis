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

//Spline
float S(float x[], float y[], float m[], float h, int j, float _x){
    float part1, part2, part3, part4;
    part1 = m[j] * (x[j+1] - _x) * (x[j+1] - _x) * (x[j+1] - _x) / (6 * h);
    part2 = m[j + 1] * (_x - x[j]) * (_x - x[j]) * (_x - x[j]) / (6 * h);
    part3 = (y[j] - m[j] * h * h / 6) * (x[j + 1] - _x) / h;
    part4 = (y[j + 1] - m[j+1] * h * h / 6) * (_x - x[j]) / h;
    return part1 + part2 + part3 + part4;
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

    //Spline
    float micro = 0.5, lambda = 0.5;
    float d[N+1];
    float _f0, _fn;
    _f0 = -32 * x[0] / ((1 + 16 * x[0] * x[0]) * (1 + 16 * x[0] * x[0]));
    _fn = -32 * x[N] / ((1 + 16 * x[N] * x[N]) * (1 + 16 * x[N] * x[N]));
    for(int i = 1; i < N; i++){
        d[i] = 3 * (y[i - 1] + y[i + 1] - 2 * y[i]) / (h * h);
    }
    d[0] = 6 / h * (((y[1] - y[0]) / h) - _f0);
    d[N] = 6 / h *(_fn - (y[N] - y[N - 1]) / h);

    float a, b, c, p[N + 1], g[N + 1], r[N], m[N + 1];
    a = 2; b = c = 0.5;
    p[0] = a; g[0] = d[0] / p[0];
    for(int i = 1; i <= N; i++){
        r[i - 1] = c / p[i-1];
        p[i] = a - b * r[i-1];
        g[i] = (d[i] - b * g[i-1]) / p[i];
    }
    m[N] = g[N];
    for(int i = N-1; i >= 0; i--){
        m[i] = g[i] - r[i] * m[i + 1];
    }

    ofstream out_x("Spline-x.txt");
	ofstream out_y("Spline-y.txt");
	if(out_x.is_open() && out_y.is_open()){
		int num = 1000;
		float len = 10.0 / num; float cur_x = -5.0 - len;
        int j = 0;
		for(int i = 0; i <= num; i++){
			cur_x += len; 
            if(cur_x >= x[j+1])
                j++;
			out_x << cur_x << endl;
			out_y << S(x, y, m, h, j, cur_x) << endl;
		}
	}
	out_x.close();
	out_y.close();


	return 0;
}
