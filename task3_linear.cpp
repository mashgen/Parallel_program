#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <ctime>
#include <fstream>
#define EPS pow(10, -6)
#define A1 0
#define B1 0
#define A2 4
#define B2 3
using namespace std;

double func_u(double x, double y){
    return sqrt(4+x*y);
}

double func_q(double x, double y){
	return (x+y);
}

double func_F(double x, double y){
    return 0.25*(x*x + y*y)/(pow(4+x*y, 3/2)) + func_q(x,y)*func_u(x,y);
}

double psi_left(double x, double y){
    return -(y / 4) + 2;
}

double psi_bottom(double x, double y){
    return -(x / 4) + 2;
}

double psi_right(double x, double y){
    return y / (4*sqrt(1+y)) + 2 * sqrt(1+y) ;
}

double psi_top(double x, double y){
    return x / (2*(sqrt(4+3*x))) + sqrt(4+3*x) ;
}


double dot(double* u, double* v, double h1, double h2, int M, int N){
    double r1, r2, res = 0;
    for (int i = 0; i <= M; ++i){
        for (int j = 0; j <= N; ++j){
            r1 = (i == 0 || i == M) ? 0.5 : 1;
            r2 = (j == 0 || j == N) ? 0.5 : 1;
            if (i != 0 && i != M && j != 0 && j != N) // скалярное произведение не определено в угловых точках
                res += u[i*(N+1) + j]*v[i*(N+1) + j]*r1*r2;
        }
    }
    res *= h1*h2;
    return res;
}

double norm(double* u, double h1, double h2, int M, int N){
    return sqrt(dot(u, u, h1, h2, M, N));
}

void diff(double* res, double* u, double* v, int M, int N){       // Разность двух векторов

    for(int i = 0; i < M+1; i++)
        for (int j = 0; j < N+1; j++)
            if(i == 0 or i == M or j == 0 or j == N) // на границах совпадают
                res[i*(N+1) + j] = 0;                           
            else
                res[i*(N+1) + j] = u[i*(N+1) + j] - v[i*(N+1) + j];
}

void vector_B(double* B, int M, int N, double h1, double h2){
    for(int i = 1; i <= M-1; i++){
        for (int j = 1; j <= N-1; j++){
            // Нижняя граница j = 0
            B[i * (N + 1) + 0] = func_F(A1 + i * h1, B1) + (2 / h2) * psi_bottom(A1 + i * h1, B1);
            // Верхняя граница j = N
            B[i * (N + 1) + N] = func_F(A1 + i * h1, B2) + (2 / h2) * psi_top(A1 + i * h1, B2);
            // Правая граница i = 0
            B[0 * (N + 1) + j] = func_F(A1, B1 + j * h2) + (2 / h1) * psi_right(A1, B1 + j * h2);
            // Левая граница i = M
            B[M * (N + 1) + j] = func_F(A2, B1 + j * h2) + (2 / h1) * psi_right(A2, B1 + j * h2);
            // Внутренние точки
            B[i * (N + 1) + j] = func_F(A1 + i * h1, B1 + j * h2);
        }
    }
    double psi_00 = (h1*psi_bottom(A1, B1) + h2*psi_left(A1, B1))/(h1+h2);
    double psi_M0 = (h1*psi_bottom(A2, B1) + h2*psi_right(A2, B1))/(h1+h2);
    double psi_MN = (h1*psi_top(A2, B2) + h2*psi_right(A2, B2))/(h1+h2);
    double psi_0N = (h1*psi_top(A1, B2) + h2*psi_left(A1, B2))/(h1+h2);

    B[0 * (N + 1) + 0] = func_F(A1 + 0 * h1, B1 + 0 * h2) + (2 / h1 + 2 / h2)*psi_00;
    B[M * (N + 1) + 0] = func_F(A1 + M * h1, B1 + 0 * h2) + (2 / h1 + 2 / h2)*psi_M0;
    B[M * (N + 1) + N] = func_F(A1 + M * h1, B1 + N * h2) + (2 / h1 + 2 / h2)*psi_MN;
    B[0 * (N + 1) + N] = func_F(A1 + 0 * h1, B1 + N * h2) + (2 / h1 + 2 / h2)*psi_0N;

    return;
}

void vector_Aw(double* A, double *w, int M, int N, double h1, double h2){
    double wx, wy;
    double temp_wx, temp_wy;

    for(int i = 1; i <= M-1; i++){

        for (int j = 1; j <= N-1; j++){


            // Нижняя граница j = 0
            temp_wx = (w[(i + 1) * (N + 1) + 0] - w[i * (N + 1) + 0])/h1 - (w[i * (N + 1) + 0] - w[(i - 1) * (N + 1) + 0])/h1;
            temp_wy = (w[i * (N + 1) + 1] - w[i * (N + 1) + 0])/h2;
            A[i * (N + 1) + 0] = -(2/h2)*temp_wy - temp_wx + (func_q(A1 + i * h1, B1 + 0 * h2) + 2/h2)* w[i * (N + 1) + 0];

            // Верхняя граница j = N
            temp_wx = (w[(i + 1) * (N + 1) + N] - w[i * (N + 1) + N])/h1 - (w[i * (N + 1) + N] - w[(i - 1) * (N + 1) + N])/h1;
            temp_wy = (w[i * (N + 1) + N] - w[i * (N + 1) + (N-1)])/h2;
            A[i * (N + 1) + N] = (2/h2)*temp_wy - temp_wx + (func_q(A1 + i * h1, B1 + N * h2) + 2/h2)* w[i * (N + 1) + N];

            // Правая граница i = 0
            temp_wx = (w[1 * (N + 1) + j] - w[0 * (N + 1) + j])/h1;
            temp_wy = (w[0 * (N + 1) + j + 1] - w[0 * (N + 1) + j])/h2 - (w[0 * (N + 1) + j] - w[0 * (N + 1) + j - 1])/h2;
            A[0 * (N + 1) + j] = -(2/h1)*temp_wx - temp_wy + (func_q(A1 + 0 * h1, B1 + j * h2) + 2/h1)* w [0 * (N + 1) + j];

            // Левая граница i = M
            temp_wx = (w[M * (N + 1) + j] - w[(M-1) * (N + 1) + j])/h1;
            temp_wy = (w[M * (N + 1) + j + 1] - w[M * (N + 1) + j])/h2 - (w[M * (N + 1) + j] - w[M * (N + 1) + j - 1])/h2;            
            A[M * (N + 1) + j] = (2/h1)*temp_wx - temp_wy + (func_q(A1 + M * h1, B1 + j * h2) + 2/h1)* w [M * (N + 1) + j];

            // Внутренние точки
            wx = (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j])/h1 - (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j])/h1;
            wy = (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j])/h2 - (w[i * (N + 1) + j] - w[i * (N + 1) + j - 1])/h2;
            A[i * (N + 1) + j] = -wx/h1 - wy/h2 + func_q(A1 + i * h1, B1 + j * h2)*w [i * (N + 1) + j];
        }
    }
    temp_wx = (w[1 * (N+1) + 0] - w[0 * (N+1) + 0]) / h1;
    temp_wy = (w[0 * (N+1) + 1] - w[0 * (N+1) + 0]) / h2;
    A[0 * (N + 1) + 0] = -(2/h1)*temp_wx - (2/h2)*temp_wy + (func_q(A1 + 0 * h1, B1 + 0 * h2) + 2/h1 + 2/h2)*w [0 * (N + 1) + 0];

    int i = M, j = 0;
    temp_wx = (w[i * (N+1) + j] - w[(i-1) * (N+1) + j]) / h1;
    temp_wy = (w[i * (N+1) + j+1] - w[i * (N+1) + j]) / h2;
    A[M * (N + 1) + 0] = (2/h1)*temp_wx - (2/h2)*temp_wy + (func_q(A1 + i * h1, B1 + j * h2) + 2/h1 + 2/h2)*w [i * (N + 1) + j];

    i = M; j = N;
    temp_wx = (w[i * (N+1) + j] - w[(i-1) * (N+1) + j]) / h1;
    temp_wy = (w[i * (N+1) + j] - w[i * (N+1) + (j-1)]) / h2;
    A[M * (N + 1) + N] = (2/h1)*temp_wx - (2/h2)*temp_wy + (func_q(A1 + i * h1, B1 + j * h2) + 2/h1 + 2/h2)*w [i * (N + 1) + j];

    i = 0; j = N;
    temp_wx = (w[(i+1) * (N+1) + j] - w[i * (N+1) + j]) / h1;
    temp_wy = (w[i * (N+1) + j] - w[i * (N+1) + (j-1)]) / h2;
    A[0 * (N + 1) + N] = (2/h1)*temp_wx - (2/h2)*temp_wy + (func_q(A1 + i * h1, B1 + j * h2) + 2/h1 + 2/h2)*w [i * (N + 1) + j];

    return;
}



int main(){

    double time; const double micro = 1.0e-07;
    struct timeval start;
    struct timeval stop;
    gettimeofday(&start,NULL);

    int M = 500, N = 500;

    double h1 = 4.0/M;
    double h2 = 3.0/N;

    double *w = new double [(M + 1) * (N + 1)];

    // Начальное значение вектора w 
    for (int i = 0; i < (M + 1) * (N + 1); i++){
        w[i] = 0;
    }

    for (int i = 0; i < M+1; i++){
        w[i*(N+1) + 0] = func_u(A1 + i* h1, B1); 
        w[i*(N+1) + N] = func_u(A1+i*h1, B2);
    }

    for(int j = 0; j < N+1; j++){
        w[0*(N+1) + j] = func_u(A1, B1 + j*h2);
        w[M*(N+1) + j] = func_u(A2, B1 +j*h2);
    }


    double *Aw = new double [(M + 1) * (N + 1)];
    double *B = new double [(M + 1) * (N + 1)];
    double *r = new double [(M + 1) * (N + 1)];
    double *Ar = new double [(M + 1) * (N + 1)];
    double *w_pred = new double [(M + 1) * (N + 1)];
    double *diff_w_w_pred = new double [(M + 1) * (N + 1)];
    double tau, norm_diff;

    cout << "Start solving\n";

    for(int i = 0; i <= (M + 1) * (N + 1); i++){
        w_pred[i] = 0;
    }
    vector_B(B, M, N, h1, h2);


    int iter = 0;
    bool solve = false;

    while(!solve){
        iter++;
        if (iter > 1){
            for (int i = 0; i < (M + 1) * (N + 1); i++)
                w_pred[i] = w[i];
        }

        vector_Aw(Aw, w_pred, M, N, h1, h2);
        diff(r, Aw, B, M, N);
        vector_Aw(Ar, r, M, N, h1, h2);
        tau = dot(Ar, r, h1, h2, M, N) / dot(Ar, Ar, h1, h2, M, N);

        for (int i = 0; i <= M; i++){
            for (int j = 0; j <= N; j++){
                if(i != 0 && i != M && j != 0 && j != N)
                    w[i * (N + 1) + j] = w_pred[i * (N + 1) + j] - tau * r[i * (N+1) + j];
            }
        }
        diff(diff_w_w_pred, w, w_pred, M, N);
        //if (iter % 10 == 0)
        //    cout << iter << ' ' << norm(diff_w_w_pred, h1, h2, M, N) << endl;
        norm_diff = norm(diff_w_w_pred, h1, h2, M, N);
        if (norm_diff < EPS)
            solve = true;
        // cout << tau << endl; 
        // cout << iter << endl;
    }
    gettimeofday(&stop,NULL);

    time = (stop.tv_sec - start.tv_sec) + micro * (stop.tv_usec - start.tv_usec);
    cout << "Time = " << time << endl;

    cout << "number of iterations: " << iter << "\n";
    cout << "tau: " << tau << "\n";
    cout << "err: " << norm_diff << "\n";

    ofstream myfile;
    myfile.open("linear_result.txt");
    // myfile << M << "\n";

    // myfile << N << "\n";
    for(int i = 0; i <= M; i++){
        for(int j = 0; j <= N; j++){
            //myfile << func_u((A1 + i * h1), (B1 + j * h2)) << ' ' <<  w[i * (N + 1) + j] << "\n";
            myfile << w[i * (N + 1) + j] << "\n";
        }
        myfile << "\n";
    }
    myfile.close ();


    delete[] w;
    delete[] w_pred;
    delete[] diff_w_w_pred;
    delete[] Aw;
    delete[] B;
    delete[] r;
    delete[] Ar;
    return 0;
}