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

void print_task() {
	cout << "Variant 17:\n\n";
	cout << "8 * x1 + 8 * x2 - 5 * x3 - 8 * x4 = 13\n";
	cout << "8 * x1 - 5 * x2 + 9 * x3 - 8 * x4 = 38\n";
	cout << "5 * x1 - 4 * x2 - 6 * x3 - 2 * x4 = 14\n";
	cout << "8 * x1 + 3 * x2 + 6 * x3 + 6 * x4 = -95\n\n";
}

// Матрица U получается путём гуассовский перестановок (маскимум выносится наверх) и приведению  к верхней треугольной
// Матрица L получается из коффициентов, на которые домножаются строчки из преобразований матрицы U и такой же 
//																									перестановке строк (нижняя треуг.)
// Матрица P получается из единичной путём перестаноки строк такоим же образом, как в преобразовании матрицы U
void LU_decomp(vector<vector<double>> mtrx_A, vector<vector<double>>& mtrx_L, vector<vector<double>>& mtrx_U, vector<vector<double>>& mtrx_P) {
	int n = mtrx_A.size();
	for (int i = 0;i < n;i++) {
		mtrx_P[i][i] = 1;
	}
	mtrx_U = mtrx_A;
	double cof;
	pair<double, int> mx;
	for (int i = 0;i < n;i++) {
		mx = { mtrx_U[i][i], i };
		for (int j = i + 1;j < n;j++) {
			if (mtrx_U[j][i] > mx.first) {
				mx = { mtrx_U[j][i], j };
			}
		}
		swap(mtrx_U[mx.second], mtrx_U[i]);
		swap(mtrx_P[mx.second], mtrx_P[i]);
		swap(mtrx_L[mx.second], mtrx_L[i]);
		for (int j = i + 1;j < n;j++) {
			cof = mtrx_U[j][i] / mtrx_U[i][i];
			mtrx_L[j][i] = cof;
			for (int k = i; k < n; k++) {
				mtrx_U[j][k] = mtrx_U[j][k] - cof * mtrx_U[i][k];
			}
		}
		mtrx_L[i][i] = 1;
	}
}

// вывод матриц
void print_mtrx(vector<vector<double>> mtrx_A, vector<vector<double>>& mtrx_L, vector<vector<double>>& mtrx_U, vector<vector<double>>& mtrx_P) {
	int n = mtrx_A.size();
	cout << "Matrix A is:\n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout << mtrx_A[i][j] << ' ';
		}
		cout << '\n';
	}
	cout << '\n';
	cout << "Matrix P is:\n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout << mtrx_P[i][j] << ' ';
		}
		cout << '\n';
	}
	cout << '\n';
	cout << "Matrix L is:\n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout << mtrx_L[i][j] << ' ';
		}
		cout << '\n';
	}
	cout << '\n';
	cout << "Matrix U is:\n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout << mtrx_U[i][j] << ' ';
		}
		cout << '\n';
	}
	cout << '\n';
}

void get_k(vector<vector<double>> mtrx_P, vector<double>& k, vector<double> b) {   // выведено самостоятельно
	int n = mtrx_P.size();
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			if (mtrx_P[i][j] != 0) {
				k[i] = b[j];
			}
		}
	}
}

void get_z(vector<vector<double>> mtrx_L, vector<double>& z, vector<double> b) {   // из методички
	int n = mtrx_L.size();
	for (int i = 0;i < n;i++) {
		z[i] = b[i];
		for (int j = 0;j < i;j++) {
			z[i] = z[i] - mtrx_L[i][j] * z[j];
		}
	}
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
	cout << "Column x is:\n\n";
	for (int i = 0;i < n;i++) {
		cout << x[i] << '\n';
	}
	cout << '\n';
}

double calc_det(vector<vector<double>> mtrx_A) {    // приведение к верхней треугольной и перемножение элементов на диагонали
	return mtrx_A[0][0] * mtrx_A[1][1] * mtrx_A[2][2] * mtrx_A[3][3];
}

double calc_det_3x3(vector<vector<double>> mtrx_A) {   // правило треугольника
	double d;
	d = mtrx_A[0][0] * mtrx_A[1][1] * mtrx_A[2][2] +
		mtrx_A[1][0] * mtrx_A[2][1] * mtrx_A[0][2] +
		mtrx_A[0][1] * mtrx_A[1][2] * mtrx_A[2][0] -
		mtrx_A[0][2] * mtrx_A[1][1] * mtrx_A[2][0] -
		mtrx_A[1][2] * mtrx_A[2][1] * mtrx_A[0][0] -
		mtrx_A[1][0] * mtrx_A[0][1] * mtrx_A[2][2];
	return d;
}

void calc_reverse_mtrx(vector<vector<double>> mtrx_A) {   // нахождение матрицы С, её транспонирование и деление
	vector<vector<double>> C;							  // всех её элементов на определитель матрицы А
	double minus;
	int n = mtrx_A.size();
	double det = calc_det(mtrx_A);
	for (int i = 0;i < n;i++) {
		vector<double> tmp;
		for (int j = 0;j < n;j++) {
			vector<vector<double>> Cij;
			for (int k = 0;k < n;k++) {
				if (k == i) {
					continue;
				}
				vector<double> cnt;
				for (int z = 0;z < n;z++) {
					if (z == j) {
						continue;
					}
					cnt.push_back(mtrx_A[k][z]);
				}
				Cij.push_back(cnt);
			}
			if ((i + j) % 2 == 1) {
				minus = -1;
			}
			else {
				minus = 1;
			}
			tmp.push_back(calc_det_3x3(Cij) / det * minus);
		}
		C.push_back(tmp);
	}
	double tmp;
	for (int i = 0;i < n;i++) {
		for (int j = i + 1;j < n;j++) {
			tmp = C[i][j];
			C[i][j] = C[j][i];
			C[j][i] = tmp;
		}
	}
	cout << "Matrix A reversed is:\n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout << C[i][j] << ' ';
		}
		cout << '\n';
	}
	cout << '\n';
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	vector<vector<double>> mtrx_A(4, vector<double>(4));
	print_task();
	//cout.precision(3);
	//cout << fixed;
	mtrx_A[0][0] = 8; mtrx_A[0][1] = 8; mtrx_A[0][2] = -5; mtrx_A[0][3] = -8;
	mtrx_A[1][0] = 8; mtrx_A[1][1] = -5; mtrx_A[1][2] = 9; mtrx_A[1][3] = -8;
	mtrx_A[2][0] = 5; mtrx_A[2][1] = -4; mtrx_A[2][2] = -6; mtrx_A[2][3] = -2;
	mtrx_A[3][0] = 8; mtrx_A[3][1] = 3; mtrx_A[3][2] = 6; mtrx_A[3][3] = 6;
	vector<double> b(4);
	b[0] = 13; b[1] = 38; b[2] = 14; b[3] = -95;
	vector<vector<double>> mtrx_L(4, vector<double>(4));
	vector<vector<double>> mtrx_U(4, vector<double>(4));
	vector<vector<double>> mtrx_P(4, vector<double>(4));
	LU_decomp(mtrx_A, mtrx_L, mtrx_U, mtrx_P);
	print_mtrx(mtrx_A, mtrx_L, mtrx_U, mtrx_P);
	vector<double> k(4);
	get_k(mtrx_P, k, b);
	vector<double> z(4);
	get_z(mtrx_L, z, k);
	vector<double> x(4);
	get_x(mtrx_U, x, z);
	cout << "Determinant of matrix A is: " << calc_det(mtrx_U) << "\n\n";
	calc_reverse_mtrx(mtrx_A);
}
