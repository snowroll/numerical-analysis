#include <iostream>
using namespace std;
#define N 100

double Hilbert(double A[N][N], int n){
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            A[i][j] = 1.0 / (i + j + 1); 
}

double get_det(double A[N][N], int n){
    if(n == 1)
        return A[0][0];
    double ans = 0; 
    double temp[N][N] = {0.0};
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n - 1; j++)
            for(int k = 0; k < n - 1; k++)
                temp[j][k] = A[j + 1][(k >= i) ? k+1 : k];
        double t = get_det(temp, n - 1);
        if(i % 2 == 0)
            ans += A[0][i] * t;
        else
            ans -=  A[0][i] * t;
    }
    return ans;
}

void  get_star(double A[N][N], double ans[N][N], int n){
    if(n == 1){
        ans[0][0] = 1;
        return;
    }
    double temp[N][N];
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n-1; k++)
                for(int t = 0; t < n - 1; t++)
                    temp[k][t] = A[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
            ans[j][i]  =  get_det(temp, n - 1);
            if((i + j) % 2 == 1)
                ans[j][i] = -ans[j][i];
        }
}

bool Inverse(double src[N][N], double des[N][N], int n){
    double detA = get_det(src, n);
    double star_matrix[N][N];
    if(detA == 0)
        return false;
    get_star(src, star_matrix, n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            des[i][j] = star_matrix[i][j] / detA;
    return true;
}

double Infinite(double A[N][N], int n){
    double ans = -10000;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            ans = max(ans, A[i][j]);
    return ans;
}

double Infinite_array(double A[N], int n){
    double ans = -10000;
    for(int i = 0; i < n; i++)
        ans = max(ans, A[i]);
    return ans;
}

void ALU(double a[N][N], double b[N], double x[N], int n)
{
    double l[n][n] = { 0 };
    double u[n][n] = { 0 };
    //进行U的第一行的赋值
    for (int i = 0; i < n; i++)
        u[0][i] = a[0][i];

    //进行L的第一列的赋值
    for (int i = 1; i < n; i++)
        l[i][0] = a[i][0] / u[0][0];

    //计算U的剩下的行数和L的剩下的列数
    for (int r = 1; r < n; r++){
        for (int i = r; i < n; i++){
            double sum1 = 0;
            for (int k = 0; k < r; k++)
                sum1 += l[r][k] * u[k][i];
                //cout << "" << r << "" << sum1 << endl;
            u[r][i] = a[r][i] - sum1;
        }

        if(r != n)
        for(int i = r+1;i<n;i++){
            double sum2 = 0;
            for(int k = 0; k < r; k++)
                sum2 += l[i][k] * u[k][r];
            l[i][r] = (a[i][r] - sum2) / u[r][r];
        }
    }

    double y[n] = { 0 };
    y[0] = b[0];
    for(int i = 1; i < n; i++){
        double sum3 = 0;
        for(int k = 0; k < i; k++)
            sum3 += l[i][k] * y[k];
        y[i] = b[i] - sum3;
    }

    x[n - 1] = y[n - 1] / u[n - 1][n - 1];
    for(int i = n - 2; i >= 0; i--){
        double sum4 = 0;
        for(int k = i + 1; k < n; k++)
            sum4 += u[i][k] * x[k];
        x[i] = (y[i] - sum4) / u[i][i];
    }
    for(int i = 0; i < n; i++)
         cout << "x[" << i + 1 << "]=" << x[i] << endl;
    return;
}

int main(){
    int n = -1;
    while(n != 0){
    double Hil_matrix[N][N] = {0.0}, Inver_matrix[N][N], b[N];
    double infi_src = 0, infi_des = 0;
    cout << "input n :" << endl;
    cin >> n;
    Hilbert(Hil_matrix, n);
    
    //Q1 求条件数
    /*Inverse(Hil_matrix, Inver_matrix, n);
    infi_src = Infinite(Hil_matrix, n); infi_des = Infinite(Inver_matrix, n);
    cout << infi_src << "   " << infi_des << "   " << infi_des * infi_src << endl; 
    */
    

    //Q2
    double x[N] = {0}, distur_flag = 0;  //distur_flag 0 没有扰动 1 有10-7的扰动 
    for(int i = 0; i < n; i++){
        double ans = 0;
        for(int j = 0; j < n; j++)
            ans += Hil_matrix[i][j];
        if(distur_flag)
            b[i] = ans - 0.0000001;
        else
            b[i] = ans;
        //cout << b[i] << endl;
    }
    ALU(Hil_matrix, b, x, n);

    double r[N] = {0};
    for(int i = 0; i < n; i++){
        double ans = 0;
        for(int j = 0; j < n; j++)
            ans += Hil_matrix[i][j] * x[j];
        r[i] = ans;
    }
    for(int i = 0; i < n; i++)
        r[i] = b[i] - r[i];
    double deta_x = 0, deta_r = 0;
    deta_r = Infinite_array(r, n) - 1;
    deta_x = Infinite_array(x, n) - 1;
    
    cout << "deta_r  " << deta_r << endl;
    cout << "deta_x  " << deta_x << endl;
    
    }
    return 0;
}
