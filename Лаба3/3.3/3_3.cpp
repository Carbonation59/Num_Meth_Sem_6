#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <queue>
#include <stack>
#include <numeric>
#include <set>

using namespace std;

double X[6] = { -3, -2, -1, 0, 1, 2 };
double Y[6] = { -2.9502, -1.8647, -0.63212, 1, 3.7183, 9.3891 };

//double X[6] = { 0, 1.7, 3.4, 5.1, 6.8, 8.5 };
//double Y[6] = { 0, 1.3038, 1.8439, 2.2583, 2.6077, 2.9155 };

vector<vector<double>> multipl(vector<vector<double>> a, vector<vector<double>> b) {
	int n = a.size();
	int n1 = a[0].size();
	vector<vector<double>> c(n, vector<double>(n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			for (int k = 0;k < n1;k++) {
				c[i][j] = c[i][j] + a[i][k] * b[k][j];
			}
		}
	}
	return c;
}

vector<double> multipl1(vector<vector<double>> a, double Y[6]) {
	int n = a.size();
	int n1 = a[0].size();
	vector<double> c(n);
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n1;j++) {
			c[i] = c[i] + a[i][j] * Y[j];
		}
	}
	return c;
}

vector<vector<double>> transpos(vector<vector<double>> a) {
	int n = a.size();
	vector<vector<double>> a1(a[0].size(), vector<double>(n));
	for (int i = 0;i < a1.size();i++) {
		for (int j = 0;j < n;j++) {
			a1[i][j] = a[j][i];
		}
	}
	return a1;
}

void get_x(vector<vector<double>> mtrx_U, vector<double>& x, vector<double> z) {   // из методички
	double cnt;
	int n = mtrx_U.size();
	for (int i = n - 1;i > -1;i--) {
		cnt = 0;
		for (int j = i + 1;j < n;j++) {
			cnt = cnt + mtrx_U[i][j] * x[j];
		}
		x[i] = (z[i] - cnt) / mtrx_U[i][i];
	}
}

vector<double> SLAY(vector<vector<double>> A, vector<double> b) {
	double mx;
	int ind;
	int sw = 0;
	int n = A.size();
	vector<vector<double>> U(n, vector<double>(n));
	vector<vector<double>> L(n, vector<double>(n));
	vector<double> w(n);
	vector<double> z(n);
	vector<double> x(n);
	vector<int> P(n);

	U = A;

	for (int i = 0; i < n; i++) {
		P[i] = i;
	}

	for (int j = 0; j < n; j++) {
		mx = abs(U[j][j]);				//поиск максимума
		ind = j;
		for (int i = j + 1; i < n; i++) {
			if (abs(U[i][j]) > mx) {
				mx = abs(U[i][j]);
				ind = i;
			}
		}

		if (j != ind) {					//замена
			sw++;
			swap(P[j], P[ind]);
			for (int i = j; i < n; i++) {
				swap(U[j][i], U[ind][i]);
			}
			for (int i = 0; i < j; i++) {
				swap(L[j][i], L[ind][i]);
			}
		}

		L[j][j] = 1;					//заполняем L(нижнетреугольную матрицу)
		for (int i = j + 1; i < n; i++) {
			L[i][j] = U[i][j] / U[j][j];
		}

		for (int i = j + 1; i < n; i++) {	//изменяем U(превращаем А в верхнетреугольную)
			for (int k = j; k < n; k++) {
				//cout << i<<" "<<k<<" "<<U[i][k]<<" "<< L[i][j]<<" "<<U[j][k]<<" "<<L[i][j]*U[j][k]<<'\n';
				U[i][k] -= L[i][j] * U[j][k];
			}
		}
	}
	/*cout << '\n';
	cout << "Матрица L(нижнетреугольная):\n";
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			cout << L[i][j] << " ";
		}
		cout << '\n';
	}
	cout << '\n';
	cout << "Матрица U(верхнетреугольная):\n";
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			cout << U[i][j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";*/
	for (int i = 0; i < n; i++) {		// Находим решение
		w[i] = b[P[i]];
	}
	for (int i = 0; i < n; i++) {
		z[i] = w[i];
		for (int j = 0; j < i; j++) {
			z[i] -= z[j] * L[i][j];
		}
	}
	for (int i = n - 1; i >= 0; i--) {
		x[i] = z[i];
		for (int j = i + 1; j < n; j++) {
			x[i] -= x[j] * U[i][j];
		}
		x[i] /= U[i][i];
	}
	/*for (int i = 0; i < n; i++){
		cout << fixed << setprecision(4) << x[i] << " ";
	}
	cout << '\n';*/
	return x;
}

void solve_mnk1(int n) {
	n++;
	vector<vector<double>> P1(n, vector <double>(2));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < 2;j++) {
			P1[i][j] = pow(X[i], j);
		}
	}
	vector<vector<double>> G = multipl(transpos(P1), P1);
	vector<double> z = multipl1(transpos(P1), Y);
	vector<double> a = SLAY(G, z);
	cout.precision(5);
	cout << fixed;
	cout << "Approximating polynomial of the 1st degree is:\n\n";
	cout << "Сoefficients is: " << a[0] << ' ' << a[1] << "\n\n";
	cout << "┌───────────────┬────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐\n";
	cout << "│       i       │      0     │     1      │     2      │     3      │     4      │     5      │\n";
	cout << "│───────────────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────│\n";
	cout << "│      x_i      │  ";
	for (int i = 0;i < n;i++) {
		if (X[i] >= 0) {
			cout << '+';
		}
		cout << X[i] << "  │  ";
	}
	cout << "\n│───────────────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────│\n";
	vector<double> pol_F(n);
	cout << "│   F_1(x_i)    │  ";
	double Fi = 0;
	for (int i = 0;i < n;i++) {
		pol_F[i] = a[0] + a[1] * X[i];
		Fi = Fi + (pol_F[i] - Y[i]) * (pol_F[i] - Y[i]);
		if (pol_F[i] >= 0) {
			cout << '+';
		}
		cout << pol_F[i] << "  │  ";
	}
	cout << "\n└───────────────┴────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘\n\n";
	cout << "Sum of Squared Errors is: " << Fi << "\n\n\n";

	vector<vector<double>> P2(n, vector <double>(3));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < 3;j++) {
			P2[i][j] = pow(X[i], j);
		}
	}
	G = multipl(transpos(P2), P2);
	z = multipl1(transpos(P2), Y);
	a = SLAY(G, z);
	cout << "Approximating polynomial of the 2nd degree is:\n\n";
	cout << "Сoefficients is: " << a[0] << ' ' << a[1] << ' ' << a[2] << "\n\n";
	cout << "┌───────────────┬────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐\n";
	cout << "│       i       │      0     │     1      │     2      │     3      │     4      │     5      │\n";
	cout << "│───────────────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────│\n";
	cout << "│      x_i      │  ";
	for (int i = 0;i < n;i++) {
		if (X[i] >= 0) {
			cout << '+';
		}
		cout << X[i] << "  │  ";
	}
	cout << "\n│───────────────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────│\n";
	cout << "│   F_2(x_i)    │  ";
	Fi = 0;
	for (int i = 0;i < n;i++) {
		pol_F[i] = a[0] + a[1] * X[i] + a[2] * X[i] * X[i];
		Fi = Fi + (pol_F[i] - Y[i]) * (pol_F[i] - Y[i]);
		if (pol_F[i] >= 0) {
			cout << '+';
		}
		cout << pol_F[i] << "  │  ";
	}
	cout << "\n└───────────────┴────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘\n\n";
	cout << "Sum of Squared Errors is: " << Fi << "\n\n\n";
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	solve_mnk1(5);
}
