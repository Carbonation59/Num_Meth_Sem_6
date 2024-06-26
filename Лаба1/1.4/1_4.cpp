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
#include <fstream>

using namespace std;

vector<vector<double>> multipl(vector<vector<double>> a, vector<vector<double>> b) {
	int n = a.size();
	vector<vector<double>> c(n, vector<double> (n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			for (int k = 0;k < n;k++) {
				c[i][j] = c[i][j] + a[i][k] * b[k][j];
			}
		}
	}
	return c;
}

vector<vector<double>> transpos(vector<vector<double>> a) {
	int n = a.size();
	double k;
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < i;j++) {
			k = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = k;
		}
	}
	return a;
}

double rotate_method(vector<vector<double>> &A, vector<vector<double>> &U_res) {
	int n = A.size();
	double mx = 0;
	pair<int, int> mx_i_j;
	for (int i = 0; i < n;i++) {
		for (int j = i + 1; j < n;j++) {
			if (abs(A[i][j]) > abs(mx)) {
				mx = A[i][j];
				mx_i_j.first = i;
				mx_i_j.second = j;
			}
		}
	}
	double phi;
	if (A[mx_i_j.first][mx_i_j.first] == A[mx_i_j.second][mx_i_j.second]) {
		phi = acos(-1) / 4;
	}
	else {
		phi = atan((2 * mx) / (A[mx_i_j.first][mx_i_j.first] - A[mx_i_j.second][mx_i_j.second])) / 2;
	}
	double cos_phi = cos(phi);
	double sin_phi = sin(phi);
	vector<vector<double>> U_cur(n, vector<double> (n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			if ((i == mx_i_j.first && j == mx_i_j.first) || (i == mx_i_j.second && j == mx_i_j.second)) {
				U_cur[i][j] = cos_phi;
			} 
			else if (i == mx_i_j.first && j == mx_i_j.second) {
				U_cur[i][j] = -sin_phi;
			}
			else if (i == mx_i_j.second && j == mx_i_j.first) {
				U_cur[i][j] = sin_phi;
			}
			else if (i == j) {
				U_cur[i][j] = 1;
			}
		}
	}
	A = multipl(transpos(U_cur), A);
	A = multipl(A, U_cur);
	double sum = 0;
	for (int i = 0; i < n;i++) {
		for (int j = i + 1; j < n;j++) {
			sum = sum + A[i][j] * A[i][j];
		}
	}
	U_res = multipl(U_res, U_cur);
	return sqrt(sum);
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	ifstream in("C:\\Users\\Artemizer\\Desktop\\6 сем. Числ. методы\\Лаба1\\1.4\\test.txt");
	int n;
	double eps;
	in >> n >> eps;
	vector<vector<double>> A(n, vector<double>(n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			in >> A[i][j];
		}
	}
	vector<vector<double>> U_res(n, vector<double>(n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			if (i == j) {
				U_res[i][j] = 1;
			}
		}
	}
	int k = 1;
	double g = rotate_method(A, U_res);
	while (g > eps) {
		g = rotate_method(A, U_res);
		k++;
	}
	cout << "matrix A:\n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout << A[i][j] << "   ";
		}
		cout << '\n';
	}
	cout << "\nlambdas is: ";
	for (int i = 0;i < n;i++) {
		cout << A[i][i] << ' ';
	}
	cout << "\n\nnumber of x is: " << n;
	cout << "\n\nx is: \n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout << U_res[i][j] << "   ";
		}
		cout << '\n';
	}
	cout << "\n\ncheck: ";
	double cnt;
	for (int i = 0;i < n;i++) {
		for (int j = i + 1;j < n;j++) {
			cnt = 0;
			for (int k = 0;k < n;k++) {
				cnt = cnt + U_res[k][i] * U_res[k][j];
			}
			cout << cnt << ' ';
		}
	}
	cout << "\nnumber of iterations is: " << k << '\n';
}
