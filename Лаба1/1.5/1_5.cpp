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

int iter = 0;

vector<vector<double>> multipl(vector<vector<double>> a, vector<vector<double>> b) {
	int n = a.size();
	vector<vector<double>> c(n, vector<double>(n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			for (int k = 0;k < n;k++) {
				c[i][j] = c[i][j] + a[i][k] * b[k][j];
			}
		}
	}
	return c;
}

double get_sign(double a) {
	if (a > 0) {
		return 1;
	}
	return -1;
}

vector<vector<double>> get_mtrx(vector<double> v) {
	int n = v.size();
	vector<vector<double>> mtrx_vv(n, vector<double>(n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			mtrx_vv[i][j] = v[i] * v[j];
		}
	}
	return mtrx_vv;
}

double get_coef(vector<double> v) {
	int n = v.size();
	double coef = 0;
	for (int i = 0;i < n;i++) {
		coef = coef + v[i] * v[i];
	}
	return coef;
}

void QR_decomp(vector<vector<double>>& A, vector<vector<double>>& Q, int k) {
	int n = A.size();
	vector<double> v(n);
	for (int i = k;i < n;i++) {
		if (i == k) {
			for (int j = i;j < n;j++) {
				v[i] = v[i] + A[j][i] * A[j][i];
			}
			v[i] = A[i][i] + get_sign(A[i][i]) * sqrt(v[i]);
		}
		else {
			v[i] = A[i][k];
		}
	}
	vector<vector<double>> mtrx_vv = get_mtrx(v);
	double coef_vv = get_coef(v);
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			mtrx_vv[i][j] = (-2) * mtrx_vv[i][j] / coef_vv;
			if (i == j) {
				mtrx_vv[i][j] = 1 + mtrx_vv[i][j];
			}
		}
	}
	A = multipl(mtrx_vv, A);
	Q = multipl(Q, mtrx_vv);
}

double find_l_val(vector<vector<double>>& A, double eps, int k) {
	int n = A.size();
	double coef = 1e9;
	while (coef > eps) {
		vector<vector<double>> R = A;
		vector<vector<double>> Q(n, vector<double>(n));
		for (int i = 0;i < n;i++) {
			for (int j = 0;j < n;j++) {
				if (i == j) {
					Q[i][j] = 1;
				}
			}
		}
		for (int i = 0;i < n - 1;i++) {
			QR_decomp(R, Q, i);
		}
		A = multipl(R, Q);
		coef = 0;
		for (int i = k + 1;i < n;i++) {
			coef = coef + A[i][k] * A[i][k];
		}
		coef = sqrt(coef);
		iter++;
	}
	return A[k][k];
}

double diff(pair<double, double> a, pair<double, double> b) {
	double c1 = a.first * a.first + a.second * a.second;
	c1 = sqrt(c1);
	double c2 = b.first * b.first + b.second * b.second;
	c2 = sqrt(c2);
	return abs(c1 - c2);
}

pair<double, double> solve(vector<vector<double>>& A, int k) {
	pair<double, double> ans = { 0,0 };
	double a = 1;
	double b = A[k][k] * (-1) - A[k + 1][k + 1];
	double c = A[k][k] * A[k + 1][k + 1] - A[k][k + 1] * A[k + 1][k];
	double D = b * b - 4 * a * c;
	if (D < 0) {
		D = -D;
		ans.second = sqrt(D) / (2 * a);
	}
	ans.first = b * (-1) / (2 * a);
	return ans;
}

pair<double, double> find_ln(vector<vector<double>>& A, double eps, int k) {
	int n = A.size();
	pair<double, double> lyambda_k = { 1, 0 };
	pair<double, double> lyambda_cur = { 0, 0 };
	while (diff(lyambda_cur, lyambda_k) > eps) {
		lyambda_k = lyambda_cur;
		vector<vector<double>> R = A;
		vector<vector<double>> Q(n, vector<double>(n));
		for (int i = 0;i < n;i++) {
			for (int j = 0;j < n;j++) {
				if (i == j) {
					Q[i][j] = 1;
				}
			}
		}
		for (int i = 0;i < n - 1;i++) {
			QR_decomp(R, Q, i);
		}
		A = multipl(R, Q);
		lyambda_cur = solve(A, k);
		iter++;
	}
	return lyambda_cur;
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	ifstream in("C:\\Users\\Artemizer\\Desktop\\6 сем. Числ. методы\\Лаба1\\1.5\\test.txt");
	int n;
	double eps;
	in >> n >> eps;
	vector<vector<double>> A(n, vector<double>(n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			in >> A[i][j];
		}
	}
	cout << "Matrix A:\n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout << A[i][j] << ' ';
		}
		cout << '\n';
	}
	cout << '\n';
	vector<vector<double>> R = A;
	vector<vector<double>> Q(n, vector<double>(n));
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			if (i == j) {
				Q[i][j] = 1;
			}
		}
	}
	for (int i = 0;i < n - 1;i++) {
		QR_decomp(R, Q, i);
	}
	cout << "Matrix Q:\n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			if (abs(Q[i][j]) < 1e-15) {
				cout << 0 << ' ';
			}
			else {
				cout << Q[i][j] << ' ';
			}
		}
		cout << '\n';
	}
	cout << "\nMatrix R:\n\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			if (abs(R[i][j]) < 1e-15) {
				cout << 0 << ' ';
			}
			else {
				cout << R[i][j] << ' ';
			}
		}
		cout << '\n';
	}
	double lyambda1;
	lyambda1 = find_l_val(A, eps, 0);
	cout << "\nlyambda 1 is: " << lyambda1 << '\n';
	pair<double, double> l2_3;
	for (int i = 1;i < n;i++) {
		if (i == n - 1) {
			cout << "lyambda " << i + 1 << " is: " << find_l_val(A, eps, i) << '\n';
			break;
		}
		l2_3 = find_ln(A, eps, i);
		if (l2_3.second == 0) {
			cout << "lyambda " << i + 1 << " is: " << find_l_val(A, eps, i) << '\n';
			//cout << "lyambda " << i + 2 << " is: " << find_l_val(A, eps, i + 1) << '\n';
		}
		else {
			cout << "lyambda " << i + 1 << " is: " << l2_3.first << ' ' << -l2_3.second << '\n';
			cout << "lyambda " << i + 2 << " is: " << l2_3.first << ' ' << l2_3.second << '\n';
			i++;
		}
	}
	cout << "\nnumber of iterations is: " << iter << '\n';
}
