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
	cout << "-19 * x1 + 2 * x2 - x3 - 8 * x4 = 38\n";
	cout << "2 * x1 + 14 * x2 - 4 * x4 = 20\n";
	cout << "6 * x1 - 5 * x2 - 20 * x3 - 6 * x4 = 52\n";
	cout << "-6 * x1 + 4 * x2 - 2 * x3 + 15 * x4 = 43\n\n";
}

const double eps = 1e-8;

vector<double> comp(vector<vector<double>> A, vector<double> a) {
	int n = a.size();
	vector<double> c(n);
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			c[i] = c[i] + A[i][j] * a[j];
		}
	}
	return c;
}

vector<double> sum(vector<double> a, vector<double> b) {
	int n = a.size();
	vector<double> c(n);
	for (int i = 0;i < n;i++) {
		c[i] = a[i] + b[i];
	}
	return c;
}

double diff(vector<double> a, vector<double> b) {
	double mx = -1;
	int n = a.size();
	for (int i = 0;i < n;i++) {
		mx = max(abs(a[i] - b[i]), mx);
	}
	return mx;
}

double mtrx_norm(vector<vector<double>> mtrx_A) {
	double sum = 1;
	int n = mtrx_A.size();
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			sum = sum + mtrx_A[i][j] * mtrx_A[i][j];
		}
	}
	sum = sqrt(sum);
	return sum;
}

int simp_iter_method(vector<vector<double>> mtrx_A, vector<double> b, vector<double>& x) {
	int n = mtrx_A.size();
	double alpha = mtrx_norm(mtrx_A);
	double coef;
	if (alpha > 1) {
		coef = alpha / (1 - alpha);
	}
	else {
		coef = 1;
	}
	coef = abs(coef);
	for (int i = 0;i < n;i++) {
		b[i] = b[i] / mtrx_A[i][i];
		for (int j = 0;j < n;j++) {
			if (i == j) {
				continue;
			}
			else {
				mtrx_A[i][j] = mtrx_A[i][j] * (-1) / mtrx_A[i][i];
			}
		}
		mtrx_A[i][i] = 0;
	}
	vector<double> x0 = b;
	x = sum(b, comp(mtrx_A, x0));
	double mx =  coef * diff(x0, x);
	int k = 1;
	while (mx > eps) {
		x0 = x;
		x = sum(b, comp(mtrx_A, x0));
		mx = coef * diff(x0, x);
		k++;
	}
	return k;
}

int zeidel_method(vector<vector<double>> mtrx_A, vector<double> b, vector<double>& x) {
	int n = mtrx_A.size();
	double alpha = mtrx_norm(mtrx_A);
	double coef;
	vector<vector<double>> mtrx_C = mtrx_A;
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < i;j++) {
			mtrx_C[i][j] = 0;
		}
	}
	double alpha1 = mtrx_norm(mtrx_C);
	if (alpha > 1) {
		coef = alpha1 / (1 - alpha);
	}
	else {
		coef = 1;
	}
	coef = abs(coef);
	for (int i = 0;i < n;i++) {
		b[i] = b[i] / mtrx_A[i][i];
		for (int j = 0;j < n;j++) {
			if (i == j) {
				continue;
			}
			else {
				mtrx_A[i][j] = mtrx_A[i][j] * (-1) / mtrx_A[i][i];
			}
		}
		mtrx_A[i][i] = 0;
	}
	vector<double> x0 = b;
	x = b;
	for (int i = 0;i < n;i++) {
		double cur = 0;
		for (int j = 0;j < n;j++) {
			if (i == j) {
				cur = cur + b[i];
			}
			else {
				cur = cur + x[j] * mtrx_A[i][j];
			}
		}
		x[i] = cur;
	}
	int k = 1;
	double mx = coef * diff(x0, x);
	while (mx > eps) {
		x0 = x;
		for (int i = 0;i < n;i++) {
			double cur = 0;
			for (int j = 0;j < n;j++) {
				if (i == j) {
					cur = cur + b[i];
				}
				else {
					cur = cur + x[j] * mtrx_A[i][j];
				}
			}
			x[i] = cur;
		}
		mx = coef * diff(x0, x);
		k++;
	}
	return k;
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	vector<vector<double>> mtrx_A(4, vector<double>(4));
	print_task();
	mtrx_A[0][0] = -19; mtrx_A[0][1] = 2; mtrx_A[0][2] = -1; mtrx_A[0][3] = -8;
	mtrx_A[1][0] = 2; mtrx_A[1][1] = 14; mtrx_A[1][2] = 0; mtrx_A[1][3] = -8;
	mtrx_A[2][0] = 6; mtrx_A[2][1] = -5; mtrx_A[2][2] = -20; mtrx_A[2][3] = -6;
	mtrx_A[3][0] = -6; mtrx_A[3][1] = 4; mtrx_A[3][2] = -2; mtrx_A[3][3] = 15;
	vector<double> b(4);
	b[0] = 38; b[1] = 20; b[2] = 52; b[3] = 43;
	vector<double> x(4);
	int k = simp_iter_method(mtrx_A, b, x);
	cout << "Simple iteration method:\n\n";
	cout << "Number of iteration is: " << k << "\n\n";
	cout << "Column x is:\n\n";
	for (int i = 0;i < 4;i++) {
		cout << x[i] << '\n';
	}
	cout << "\n\n";
	k = zeidel_method(mtrx_A, b, x);
	cout << "Zeidel method:\n\n";
	cout << "Number of iteration is: " << k << "\n\n";
	cout << "Column x is:\n\n";
	for (int i = 0;i < 4;i++) {
		cout << x[i] << '\n';
	}
	cout << '\n';
}
