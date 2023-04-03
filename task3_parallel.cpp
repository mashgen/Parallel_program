#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <fstream>

#include <mpi.h>
using namespace std;

#define EPS 5*pow(10, -6)


double func_u(double x, double y) {
    return sqrt(4 + x * y);
}

double func_q(double x, double y) {
    return (x + y);
}

double func_F(double x, double y) {
    return 0.25 * (x * x + y * y) / (pow(4 + x * y, 3 / 2)) + func_q(x, y) * func_u(x, y);
}

double psi_left(double x, double y) {
    return -(y / 4) + 2;
}

double psi_bottom(double x, double y) {
    return -(x / 4) + 2;
}

double psi_right(double x, double y) {
    return y / (4 * sqrt(1 + y)) + 2 * sqrt(1 + y);
}

double psi_top(double x, double y) {
    return x / (2 * (sqrt(4 + 3 * x))) + sqrt(4 + 3 * x);
}


double dot(double* u, double* v, double h1, double h2, int M, int N) {
    double r1, r2, res = 0;
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            r1 = (i == 0 || i == M) ? 0.5 : 1;
            r2 = (j == 0 || j == N) ? 0.5 : 1;
            if (i != 0 && i != M && j != 0 && j != N) // скалярное произведение не определено в угловых точках
                res += u[i * (N + 1) + j] * v[i * (N + 1) + j] * r1 * r2;
        }
    }
    res *= h1 * h2;
    return res;
}

double norm(double* u, double h1, double h2, int M, int N) {
    return sqrt(dot(u, u, h1, h2, M, N));
}

void diff(double* res, double* u, double* v, int M, int N) {       // Разность двух векторов

    for (int i = 0; i < M + 1; i++)
        for (int j = 0; j < N + 1; j++)
            if (i == 0 or i == M or j == 0 or j == N) // на границах совпадают
                res[i * (N + 1) + j] = 0;
            else
                res[i * (N + 1) + j] = u[i * (N + 1) + j] - v[i * (N + 1) + j];
}

void vector_B(double* B, int M, int N, double h1, double h2, double x_0, double y_0) {
    double A1 = x_0, B1 = y_0;
    for (int i = 1; i <= M - 1; i++) {
        for (int j = 1; j <= N - 1; j++) {
            // Нижняя граница j = 0
            B[i * (N + 1) + 0] = func_F(A1 + i * h1, B1) + (2 / h2) * psi_bottom(A1 + i * h1, B1);
            // Верхняя граница j = N
            B[i * (N + 1) + N] = func_F(A1 + i * h1, B1 + j * N) + (2 / h2) * psi_top(A1 + i * h1, B1 + j * N);
            // Правая граница i = 0
            B[0 * (N + 1) + j] = func_F(A1, B1 + j * h2) + (2 / h1) * psi_right(A1, B1 + j * h2);
            // Левая граница i = M
            B[M * (N + 1) + j] = func_F(A1 + M * h1, B1 + j * h2) + (2 / h1) * psi_right(A1 + M * h1, B1 + j * h2);
            // Внутренние точки
            B[i * (N + 1) + j] = func_F(A1 + i * h1, B1 + j * h2);
        }
    }
    double psi_00 = (h1 * psi_bottom(A1, B1) + h2 * psi_left(A1, B1)) / (h1 + h2);
    double psi_M0 = (h1 * psi_bottom(A1 + M * h1, B1) + h2 * psi_right(A1 + M * h1, B1)) / (h1 + h2);
    double psi_MN = (h1 * psi_top(A1 + M * h1, B1 + N * h2) + h2 * psi_right(A1 + M * h1, B1 + N * h2)) / (h1 + h2);
    double psi_0N = (h1 * psi_top(A1, B1 + N * h2) + h2 * psi_left(A1, B1 + N * h2)) / (h1 + h2);

    B[0 * (N + 1) + 0] = func_F(A1 + 0 * h1, B1 + 0 * h2) + (2 / h1 + 2 / h2) * psi_00;
    B[M * (N + 1) + 0] = func_F(A1 + M * h1, B1 + 0 * h2) + (2 / h1 + 2 / h2) * psi_M0;
    B[M * (N + 1) + N] = func_F(A1 + M * h1, B1 + N * h2) + (2 / h1 + 2 / h2) * psi_MN;
    B[0 * (N + 1) + N] = func_F(A1 + 0 * h1, B1 + N * h2) + (2 / h1 + 2 / h2) * psi_0N;

    return;
}

void vector_Aw(double* A, double* w, int M, int N, double h1, double h2, double x_0, double y_0) {
    double A1 = x_0, B1 = y_0;
    double wx, wy;
    double temp_wx, temp_wy;

    for (int i = 1; i <= M - 1; i++) {

        for (int j = 1; j <= N - 1; j++) {


            // Нижняя граница j = 0
            temp_wx = (w[(i + 1) * (N + 1) + 0] - w[i * (N + 1) + 0]) / h1 - (w[i * (N + 1) + 0] - w[(i - 1) * (N + 1) + 0]) / h1;
            temp_wy = (w[i * (N + 1) + 1] - w[i * (N + 1) + 0]) / h2;
            A[i * (N + 1) + 0] = -(2 / h2) * temp_wy - temp_wx + (func_q(A1 + i * h1, B1 + 0 * h2) + 2 / h2) * w[i * (N + 1) + 0];

            // Верхняя граница j = N
            temp_wx = (w[(i + 1) * (N + 1) + N] - w[i * (N + 1) + N]) / h1 - (w[i * (N + 1) + N] - w[(i - 1) * (N + 1) + N]) / h1;
            temp_wy = (w[i * (N + 1) + N] - w[i * (N + 1) + (N - 1)]) / h2;
            A[i * (N + 1) + N] = (2 / h2) * temp_wy - temp_wx + (func_q(A1 + i * h1, B1 + N * h2) + 2 / h2) * w[i * (N + 1) + N];

            // Правая граница i = 0
            temp_wx = (w[1 * (N + 1) + j] - w[0 * (N + 1) + j]) / h1;
            temp_wy = (w[0 * (N + 1) + j + 1] - w[0 * (N + 1) + j]) / h2 - (w[0 * (N + 1) + j] - w[0 * (N + 1) + j - 1]) / h2;
            A[0 * (N + 1) + j] = -(2 / h1) * temp_wx - temp_wy + (func_q(A1 + 0 * h1, B1 + j * h2) + 2 / h1) * w[0 * (N + 1) + j];

            // Левая граница i = M
            temp_wx = (w[M * (N + 1) + j] - w[(M - 1) * (N + 1) + j]) / h1;
            temp_wy = (w[M * (N + 1) + j + 1] - w[M * (N + 1) + j]) / h2 - (w[M * (N + 1) + j] - w[M * (N + 1) + j - 1]) / h2;
            A[M * (N + 1) + j] = (2 / h1) * temp_wx - temp_wy + (func_q(A1 + M * h1, B1 + j * h2) + 2 / h1) * w[M * (N + 1) + j];

            // Внутренние точки
            wx = (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j]) / h1 - (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j]) / h1;
            wy = (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j]) / h2 - (w[i * (N + 1) + j] - w[i * (N + 1) + j - 1]) / h2;
            A[i * (N + 1) + j] = -wx / h1 - wy / h2 + func_q(A1 + i * h1, B1 + j * h2) * w[i * (N + 1) + j];
        }
    }
    temp_wx = (w[1 * (N + 1) + 0] - w[0 * (N + 1) + 0]) / h1;
    temp_wy = (w[0 * (N + 1) + 1] - w[0 * (N + 1) + 0]) / h2;
    A[0 * (N + 1) + 0] = -(2 / h1) * temp_wx - (2 / h2) * temp_wy + (func_q(A1 + 0 * h1, B1 + 0 * h2) + 2 / h1 + 2 / h2) * w[0 * (N + 1) + 0];

    int i = M, j = 0;
    temp_wx = (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j]) / h1;
    temp_wy = (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j]) / h2;
    A[M * (N + 1) + 0] = (2 / h1) * temp_wx - (2 / h2) * temp_wy + (func_q(A1 + i * h1, B1 + j * h2) + 2 / h1 + 2 / h2) * w[i * (N + 1) + j];

    i = M; j = N;
    temp_wx = (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j]) / h1;
    temp_wy = (w[i * (N + 1) + j] - w[i * (N + 1) + (j - 1)]) / h2;
    A[M * (N + 1) + N] = (2 / h1) * temp_wx - (2 / h2) * temp_wy + (func_q(A1 + i * h1, B1 + j * h2) + 2 / h1 + 2 / h2) * w[i * (N + 1) + j];

    i = 0; j = N;
    temp_wx = (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j]) / h1;
    temp_wy = (w[i * (N + 1) + j] - w[i * (N + 1) + (j - 1)]) / h2;
    A[0 * (N + 1) + N] = (2 / h1) * temp_wx - (2 / h2) * temp_wy + (func_q(A1 + i * h1, B1 + j * h2) + 2 / h1 + 2 / h2) * w[i * (N + 1) + j];

    return;
}

void get_domens(int* dim, int M, int N, int numproc) {
    // m по x, n по y
    int n = 0, m = 0;
    int A = M, B = N;
    while (pow(2, m + n) != numproc) {
        if (A > B) {
            m++;
            A = A / 2;
        }
        else {
            n++;
            B = B / 2;
        }
    }
    dim[0] = pow(2, m);
    dim[1] = pow(2, n);

    return;
}


int main(int argc, char** argv)
{
    int A1 = 0, B1 = 0, A2 = 4, B2 = 3;
    int M = 500, N = 500;

    double h1 = 4.0 / M;
    double h2 = 3.0 / N;

    int rank, numprocs;
    int dim[2];
    dim[0] = 0;
    dim[1] = 0;
    int ndims = 2; // двумерное разбиение

    //cout << "Start parallel" << endl;

    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm MPI_COMM_CART;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    //cout <<"rank: " << rank << endl;

    //MPI_Dims_create(numprocs, 2, dim);
    get_domens(dim, M + 1, N + 1, numprocs);

    int def_amount_x = (M + 1) / dim[0];
    int def_amount_y = (N + 1) / dim[1];
    int rem_amount_x = (M + 1) % dim[0];
    int rem_amount_y = (N + 1) % dim[1];

    int amount_x[8][8];
    int amount_y[8][8];


    // получаем все точки разбиения по x

    for (int i = 0; i < dim[0]; i++) { // идем по столбцам
        if (rem_amount_x > 0) {
            for (int j = 0; j < dim[1]; j++)
                amount_x[i][j] = def_amount_x + 1;
            rem_amount_x--;
        }
        else {
            for (int j = 0; j < dim[1]; j++)
                amount_x[i][j] = def_amount_x;
        }
    }

    // получаем все точки разбиения по y

    for (int j = 0; j < dim[1]; j++) { // идем по столбцам
        if (rem_amount_y > 0) {
            for (int i = 0; i < dim[0]; i++)
                amount_y[i][j] = def_amount_y + 1;
            rem_amount_y--;
        }
        else {
            for (int i = 0; i < dim[0]; i++)
                amount_y[i][j] = def_amount_y;
        }
    }


    //if (rank == 0)
     //   cout << rank << " : " << " DIMS " << dim[0] << " " << dim[1] << "\n";

    // нет периодичности у решетки
    int periods[2];
    periods[0] = 0;
    periods[1] = 0;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dim, periods, true, &MPI_COMM_CART);

    // для получение координат процессов
    int coords[2];
    MPI_Cart_coords(MPI_COMM_CART, rank, ndims, coords);

    int x_size, y_size;
    x_size = amount_x[coords[0]][coords[1]];
    y_size = amount_y[coords[0]][coords[1]];
    //cout << rank << " : " << " Coords " << amount_x[coords[0]][coords[1]] << "*" << amount_y[coords[0]][coords[1]] << endl;


    double x_0, y_0;
    x_0 = A1 + x_size * coords[0];
    y_0 = B1 + y_size * coords[1];

    double* w = new double[(x_size + 2) * (y_size + 2)];

    for (int i = 0; i <= x_size + 1; i++) {
        for (int j = 0; j <= y_size + 1; j++) {
            w[i * (y_size + 2) + j] = 0;
        }
    }

    int source_less_x, source_greater_x, dest_less_x, dest_greater_x;
    int source_less_y, source_greater_y, dest_less_y, dest_greater_y;
    if (coords[0] == 0) { //слева, не отправляют влево ничего и не получают оттуда
        source_less_x = dest_less_x = MPI_PROC_NULL;
        for (int j = 0; j <= y_size + 1; j++)
            w[0 * (y_size + 2) + j] = func_u(A1, B1 + j * h2);
    }
    else
        source_less_x = dest_less_x = coords[0] - 1;

    if (coords[0] == dim[0] - 1) { //справа, не отправляют  ничего и не получают оттуда
        source_greater_x = dest_greater_x = MPI_PROC_NULL;
        for (int j = 0; j <= y_size + 1; j++) {
            w[x_size * (y_size + 2) + j] = func_u(A2, B1 + j * h2);
        }
    }
    else
        source_greater_x = dest_greater_x = coords[0] + 1;

    if (coords[1] == 0) {//снизу, не отправляют  ничего и не получают оттуда
        source_less_y = dest_less_y = MPI_PROC_NULL;
        for (int i = 0; i <= x_size + 1; i++)
            w[i * (y_size + 2) + 0] = func_u(A1 + i * h1, B1);
    }
    else
        source_less_y = dest_less_y = coords[1] - 1;

    if (coords[1] == dim[1] - 1) { //сверху, не отправляют  ничего и не получают оттуда
        source_greater_y = dest_greater_y = MPI_PROC_NULL;
        for (int i = 0; i <= x_size + 1; i++)
            w[i * (y_size + 2) + y_size] = func_u(A1 + i * h1, B1);

    }
    else
        source_greater_y = dest_greater_y = coords[1] + 1;

    // получим rank получателя и от кого приходит сообщение

    int recv_less_x, send_less_x, recv_greater_x, send_greater_x;
    int recv_less_y, send_less_y, recv_greater_y, send_greater_y;

    int sendrecv_coord[2];
    sendrecv_coord[0] = 0;
    sendrecv_coord[1] = 0;

    //по координатам получаем rank, 4 раза, если есть кому присылать и от кого принимать
    if (source_less_x != MPI_PROC_NULL) {
        sendrecv_coord[0] = source_less_x;
        sendrecv_coord[1] = coords[1];
        MPI_Cart_rank(MPI_COMM_CART, sendrecv_coord, &recv_less_x);
        send_less_x = recv_less_x;
    }
    else {
        send_less_x = recv_less_x = MPI_PROC_NULL;
    }


    if (source_greater_x != MPI_PROC_NULL) {
        sendrecv_coord[0] = source_greater_x;
        sendrecv_coord[1] = coords[1];
        MPI_Cart_rank(MPI_COMM_CART, sendrecv_coord, &recv_greater_x);
        send_greater_x = recv_greater_x;
    }
    else {
        send_greater_x = recv_greater_x = MPI_PROC_NULL;
    }

    if (source_less_y != MPI_PROC_NULL) {
        sendrecv_coord[1] = source_less_y;
        sendrecv_coord[0] = coords[0];
        MPI_Cart_rank(MPI_COMM_CART, sendrecv_coord, &recv_less_y);
        send_less_y = recv_less_y;
    }
    else
        send_less_y = recv_less_y = MPI_PROC_NULL;

    if (source_greater_y != MPI_PROC_NULL) {
        sendrecv_coord[1] = source_greater_y;
        sendrecv_coord[0] = coords[0];
        MPI_Cart_rank(MPI_COMM_CART, sendrecv_coord, &recv_greater_y);
        send_greater_y = recv_greater_y;
    }
    else {
        send_greater_y = recv_greater_y = MPI_PROC_NULL;
    }

    //cout << "rank: " << rank << "send_to: " << send_greater_x << " " << send_less_x << " " << send_greater_y << " " << send_less_y << endl;


    //буферы для пересылки

    double* top_s, * top_r;
    top_s = new double[x_size];
    top_r = new double[x_size];
    double* bottom_s, * bottom_r;
    bottom_s = new double[x_size];
    bottom_r = new double[x_size];
    double* left_s, * left_r;
    left_s = new double[y_size];
    left_r = new double[y_size];
    double* right_s, * right_r;
    right_s = new double[y_size];
    right_r = new double[y_size];



    double* w_pred = new double[(x_size + 2) * (y_size + 2)];
    double* Aw = new double[(x_size + 2) * (y_size + 2)];
    double* B = new double[(x_size + 2) * (y_size + 2)];
    double* r = new double[(x_size + 2) * (y_size + 2)];
    double* Ar = new double[(x_size + 2) * (y_size + 2)];
    double* diff_w_w_pred = new double[(x_size + 2) * (y_size + 2)];

    double tau;
    double eps_local, eps_forall = 10;

    int iter = 0;
    bool solve = false;

    double start_time = MPI_Wtime();

    vector_B(B, x_size + 1, y_size + 1, h1, h2, A1 + x_size * coords[0] * h1, B1 + y_size * coords[1] * h2);

    for (int i = 0; i < x_size + 2; i++)
        for (int j = 0; j < y_size + 2; j++)
            w_pred[i * (y_size + 2) + j] = 0;

    while (!solve) {
        iter++;
        if (iter > 1) {
            for (int i = 0; i <= x_size + 1; i++)
                for (int j = 0; j <= y_size + 1; j++)
                    w_pred[i * (y_size + 2) + j] = w[i * (y_size + 2) + j];
        }
        // заполняем буферы
        for (int i = 1; i < x_size + 1; i++) {
            bottom_s[i - 1] = w[i * (y_size + 2) + 1];
            top_s[i - 1] = w[i * (y_size + 2) + y_size];
        }

        for (int j = 1; j < y_size + 1; j++) {
            right_s[j - 1] = w[(y_size + 2) + j];
            left_s[j - 1] = w[x_size * (y_size + 2) + j];
        }

        // отправляем буферы
        // по вертикали вверх отправляем, снизу получаем
        MPI_Sendrecv(top_s, x_size, MPI_DOUBLE, send_greater_y, 0, bottom_r, x_size, MPI_DOUBLE, recv_less_y, 0, MPI_COMM_CART, &status);
        // по вертикали вниз отправляем, сверху получаем
        MPI_Sendrecv(bottom_s, x_size, MPI_DOUBLE, send_less_y, 0, top_r, x_size, MPI_DOUBLE, recv_greater_y, 0, MPI_COMM_CART, &status);

        // по горизонтали вправо отправляем, слева получаем
        MPI_Sendrecv(right_s, y_size, MPI_DOUBLE, send_greater_x, 0, left_r, y_size, MPI_DOUBLE, recv_less_x, 0, MPI_COMM_CART, &status);
        // по горизонтали влево отправляем, справа получаем
        MPI_Sendrecv(left_s, y_size, MPI_DOUBLE, send_less_x, 0, right_r, y_size, MPI_DOUBLE, recv_greater_x, 0, MPI_COMM_CART, &status);


        // заполняем из буфера

        for (int i = 1; i < x_size + 1; i++) {
            w[i * (y_size + 2) + 1] = bottom_r[i - 1];
            w[i * (y_size + 2) + y_size] = top_r[i - 1];
        }

        for (int j = 1; j < y_size + 1; j++) {
            w[(y_size + 2) + j] = right_r[j - 1];
            w[x_size * (y_size + 2) + j] = left_r[j - 1];
        }

        vector_Aw(Aw, w_pred, x_size + 1, y_size + 1, h1, h2, A1 + x_size * coords[0] * h1, B1 + y_size * coords[1] * h2);
        diff(r, Aw, B, x_size + 1, y_size + 1);
        vector_Aw(Ar, r, x_size + 1, y_size + 1, h1, h2, A1 + x_size * coords[0] * h1, B1 + y_size * coords[1] * h2);
        tau = dot(Ar, r, h1, h2, x_size + 1, y_size + 1) / dot(Ar, Ar, h1, h2, x_size + 1, y_size + 1);

        for (int i = 0; i <= x_size + 1; i++) {
            for (int j = 0; j <= y_size + 1; j++) {
                if (i != 0 && i != x_size + 1 && j != 0 && j != y_size + 1)
                    w[i * (y_size + 2) + j] = w_pred[i * (y_size + 2) + j] - tau * r[i * (y_size + 2) + j];
            }
        }
        diff(diff_w_w_pred, w, w_pred, x_size + 1, y_size + 1);

        eps_local = norm(diff_w_w_pred, h1, h2, x_size + 1, y_size + 1);

        eps_forall = 0;

        MPI_Allreduce(&eps_local, &eps_forall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_CART);

        if ((rank == 0) && (iter % 500 == 0))
            cout <<"eps: " << eps_forall / 1 << " iter: " << iter <<endl;
        if (eps_forall < EPS)
            solve = true;
    }

    MPI_Barrier(MPI_COMM_CART);
    double end_time = MPI_Wtime();
    if (rank == 0) {
        cout << "TIME " << end_time - start_time << endl;
        cout << iter;
    }


    delete[] w;
    delete[] w_pred;
    delete[] diff_w_w_pred;
    delete[] Aw;
    delete[] B;
    delete[] r;
    delete[] Ar;

    delete[] top_s;
    delete[] top_r;
    delete[] bottom_s;
    delete[] bottom_r;
    delete[] left_s;
    delete[] left_r;
    delete[] right_s;
    delete[] right_r;

    MPI_Comm_free(&MPI_COMM_CART);
    MPI_Finalize();

    return 0;
}