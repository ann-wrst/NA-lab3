
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
const int N = 2;
const int M = 3;
struct point {
    int i;
    int j;
    double max;
};
double scalarProduct(double a[M], double b[M]) {
    double sp = 0;
    for (int i = 0; i < M; i++) {
        sp += a[i] * b[i];
    }
    return sp;
}
double MatrixNorm(double** ar, int m, int n) {
    double sum;
    double max = 0;
    for (int i = 0; i < n; i++) {
        sum = 0;
        for (int j = 0; j < m; j++) {
            sum += abs(ar[i][j]);
        }
        if (sum > max)
            max = sum;
    }
    return max;
}
void maxInColumn(int m, double **ar, point* temp) {
    temp->max = abs(ar[0][m]);
    temp->i = 0;
    temp->j = m;
    for (int j = 0; j < N - 1; j++) {
        if (abs(ar[j + 1][m]) > temp->max) {
            temp->max = abs(ar[j + 1][m]);
            temp->i = j + 1;
            temp->j = m;
        }
    }
}
void inverseGaussMethod(double **arD) {
    double** e;
    e = new double* [N];
    for (int i = 0; i < N; i++) {
        e[i] = new double[N];
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            e[i][j] = 0.0;
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            e[i][i] = 1;
        }
    }
    point Point[N];
    for (int s = 0; s < N - 1; s++) {
        maxInColumn(s, arD, &Point[s]);
        if (Point[s].i != s) {  
            for (int k = 0; k < N; k++) {
                swap(arD[Point[s].i][k], arD[s][k]);
                swap(e[Point[s].i][k], e[s][k]);
            }
            Point[s].i = s;
        }
        for (int k = s + 1; k < N; k++) {
            if (arD[k][Point[s].j] != 0) {
                double m = -(double)arD[k][Point[s].j] / arD[Point[s].i][Point[s].j];

                for (int i = 0; i < N; i++) {
                    e[k][i] += e[Point[s].i][i] * m;
                    arD[k][i] += arD[Point[s].i][i] * m;
                }
            }
        }
    }

    for (int i = 0; i < N; i++) {
        double diag = arD[i][i];
        for (int j = 0; j < N; j++) {
            if (diag != 0) {
                arD[i][j] /= diag;
                e[i][j] /= diag;
            }
        }
    }

    for (int j = N - 1; j >= 1; j--) {
        for (int i = j - 1; i >= 0; i--) {
            double m = -arD[i][j];
            for (int k = 0; k < N; k++) {
                arD[i][k] += arD[j][k] * m;
                e[i][k] += e[j][k] * m;
            }
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            arD[i][j] = e[i][j];
        }
    }
}

double eq1(double x,double y) {
    return sin(2 * x - y) - 1.2 * x - 0.4;
}
double eq2(double x, double y) {
    return 0.8 * x * x + 1.5 * y * y - 1;
}
double eq1DerivativeX(double x, double y) {
    return 2 * cos(2 * x - y) - 1.2;
}
double eq1DerivativeY(double x, double y) {
    return -cos(2 * x - y);
}
double eq2DerivativeX(double x, double y) {
    return 1.6*x;
}
double eq2DerivativeY(double x, double y) {
    return 3*y;
}
void newtonMethod(double appr[N], double xr[N], double eps = 1e-3) {
    double** A;
    A = new double* [N];
    for (int i = 0; i < N; i++) {
        A[i] = new double[N];
    }
    double fx[N];
    double z[N];
    int n = 1;
    int choice;
    cout << "Newton method" << endl;
    cout << "Enter 1, if you want to set a number of iterations\nEnter 2, if you want to get accuracy " << eps << endl;
    cin >> choice;
    if (choice == 1) {
        cout << "Enter the number of iterations\n";
        int iterNumber;
        cin >> iterNumber;
        cout << " " << "\t\tx\t\t  y" << endl;
        for (int i = 0; i < iterNumber; i++) {
            A[0][0] = eq1DerivativeX(appr[0], appr[1]);
            A[0][1] = eq1DerivativeY(appr[0], appr[1]);
            A[1][0] = eq2DerivativeX(appr[0], appr[1]);
            A[1][1] = eq2DerivativeY(appr[0], appr[1]);

            fx[0] = eq1(appr[0], appr[1]);
            fx[1] = eq2(appr[0], appr[1]);
            inverseGaussMethod(A);
            //умножение обратной А на fx
            for (int i = 0; i < N; i++) {
                z[i] = A[i][0] * fx[0] + A[i][1] * fx[1];
            }
            for (int i = 0; i < N; i++) {
                appr[i] = appr[i] - z[i];
            }
            cout << n;
            cout.width(20);
            cout << right << appr[0];
            cout.width(20);
            cout << right << appr[1] << endl;
            n++;
        }
    }
    else if (choice == 2) {
        cout << " " << "\t\tx\t\t  y" << endl;
        for (int i = 0; i < N; i++) {
            do {
                //заполнение матрицы якоби
                A[0][0] = eq1DerivativeX(appr[0], appr[1]);
                A[0][1] = eq1DerivativeY(appr[0], appr[1]);
                A[1][0] = eq2DerivativeX(appr[0], appr[1]);
                A[1][1] = eq2DerivativeY(appr[0], appr[1]);

                fx[0] = eq1(appr[0], appr[1]);
                fx[1] = eq2(appr[0], appr[1]);
                inverseGaussMethod(A);
                //умножение обратной А на fx
                for (int i = 0; i < N; i++) {
                    z[i] = A[i][0] * fx[0] + A[i][1] * fx[1];
                }
                for (int i = 0; i < N; i++) {
                    appr[i] = appr[i] - z[i];
                }
                cout << n;
                cout.width(20);
                cout << right << appr[0];
                cout.width(20);
                cout << right << appr[1] << endl;
                n++;

            } while (abs(appr[i] - xr[i]) > eps);

        }
    }
    else cout << "The input is incorrect\n";
}
void powerMethod(double** A, double xk[M], double eps=1e-3) {
    cout << endl << "Power method\n" << "Matrix:\n";
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    double matrixNorm = MatrixNorm(A, M, M);
    double e[M][M];
    double lambdaMax, lambdaMin;
    double temp[M];
    double m = 0;
    double prevM = 1;
    double B[M][M];
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++)
            if (i != j)
                e[i][j] = 0.0;
        e[i][i] = matrixNorm;

    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            B[i][j] = e[i][j] - A[i][j];
        }
    }
    while (abs(m - prevM) >= eps) {
        //умножение B на xk
        for (int i = 0; i < M; i++) {
            temp[i] = B[i][0] * xk[0] + B[i][1] * xk[1] + B[i][2] * xk[2];
        }
        prevM = m;
        m = scalarProduct(temp, xk) / scalarProduct(xk, xk);
        for (int i = 0; i < M; i++) {
            xk[i] = temp[i];
        }
    }
    lambdaMax = m;
    lambdaMin = matrixNorm - lambdaMax;
    cout << "Мiнiмальне власне значення: " << lambdaMin << endl;
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << fixed << setprecision(5);
    double appr[N] = { 1.25,0 };
    double xr[N] = {-0.941581562518972, 0.440256992296897};
    newtonMethod(appr,xr);
    double xk[M] = { 1.0,1.0,1.0 };
    double** A = new double* [M];
    for (int i = 0; i < M; i++) {
        A[i] = new double[M];
    }
    A[0][0] = 2.0;
    A[0][1] = 1.0;
    A[0][2] = 0.0;
    A[1][0] = 1.0;
    A[1][1] = 2.0;
    A[1][2] = 1.0;
    A[2][0] = 0.0;
    A[2][1] = 1.0;
    A[2][2] = 2.0;
    powerMethod(A, xk);
}
