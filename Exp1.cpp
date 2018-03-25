#include <iostream>
#include <math.h>
using namespace std;

int main(){
	double ln2 = 0.693147190546;
	float epsilon = 0.00005;
	//float xn = 1;
	double xn = 1;
	int k = 1;  //实际指数的k-1 为奇偶之分
	while(fabs(xn - ln2) >= epsilon){
		k++;
		if(k % 2 == 0)  //k为偶
			xn += (-1 * 1.0) / k;
		else
			xn += 1.0 / k;
	}
	cout << " n is " << k;

	return 0;
}