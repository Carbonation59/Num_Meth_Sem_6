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

double X_1[4] = { -2, -1, 0, 1 };
double X_2[4] = { -2, -1, 0.2, 1 };
double X_st = -0.5;

double func(double x) {
	return exp(x) + x;
}

double omega(double x, double X[4], int j) {
	double ans = 1;
	for (int i = 0;i < 4;i++) {
		if (j == i) {
			continue;
		}
		ans = ans * (x - X[i]);
	}
	return ans;
}

void interp_lagr(double X[4], int n) {
	cout << "Lagrange:\n\n";
	cout << "┌───────────────┬─────────────┬────────────┬────────────┬────────────────────┬────────────┐\n";
	cout << "│      i        │     x_i     │     f_i    │  w'_4(x_i) │   f_i / w'_4(x_i)  │  X* - x_i  │\n";
	cout.precision(5);
	cout << fixed;
	vector<double> v;
	for (int i = 0;i < n;i++) {
		cout << "│───────────────┼─────────────┼────────────┼────────────┼────────────────────┼────────────│\n";
		cout << "│      " << i << "        │   ";
		if (X[i] >= 0) {
			cout << '+';
		}
		cout << X[i] << "  │  ";
		if (func(X[i]) >= 0) {
			cout << '+';
		}
		cout << func(X[i]) << "  │  ";
		if (omega(X[i], X, i) >= 0) {
			cout << '+';
		}
		cout << omega(X[i], X, i) << "  │      ";
		if (func(X[i]) / omega(X[i], X, i) >= 0) {
			cout << '+';
		}
		cout << func(X[i]) / omega(X[i], X, i) << "      │  ";
		if (X_st - X[i] >= 0) {
			cout << '+';
		}
		cout << X_st - X[i] << "  │\n";
		v.push_back(func(X[i]) / omega(X[i], X, i));
	}
	cout << "└───────────────┴─────────────┴────────────┴────────────┴────────────────────┴────────────┘\n\n";
	double L = 0;
	for (int i = 0;i < n;i++) {
		L = L + v[i] * omega(X_st, X, i);
	}
	cout << "L_3(X*) = " << L << '\n';;
	cout << "y(X*) = " << func(X_st) << '\n';
	cout << "Δ(L_3(X*)) = " << abs(L - func(X_st)) << "\n\n";
}

double newt_calc(double X[4], int j, double x) {
	double ans = 1;
	for (int i = 0;i <= j;i++) {
		ans = ans * (x - X[i]);
	}
	return ans;
}

void interp_newt(double X[4], int n) {
	cout << "Newton:\n\n";
	cout << "┌───────────────┬─────────────┬────────────┐\n";
	cout << "│      i        │     x_i     │     f_i    │\n";
	cout.precision(5);
	cout << fixed;
	vector<double> v;
	for (int i = 0;i < n;i++) {
		cout << "│───────────────┼─────────────┼────────────│\n";
		cout << "│      " << i << "        │   ";
		if (X[i] >= 0) {
			cout << '+';
		}
		cout << X[i] << "  │  ";
		if (func(X[i]) >= 0) {
			cout << '+';
		}
		cout << func(X[i]) << "  │\n";
		v.push_back(func(X[i]));
	}
	cout << "└───────────────┴─────────────┴────────────┘\n\n";
	int N = v.size();
	int step;
	double a;
	vector<double> coefs;
	while (v.size() != 1) {
		cout << "┌────────────┐\n";
		cout << "│     f_i    │\n";
		vector<double> v1;
		step = N - v.size() + 1;
		for (int i = 0;i < v.size() - 1;i++) {
			cout << "│────────────│\n";
			a = (v[i] - v[i + 1]) / (X[i] - X[i + step]);
			cout << "│   " << a << "  │\n";
			v1.push_back(a);
		}
		v = v1;
		coefs.push_back(v[0]);
		cout << "└────────────┘\n\n";
	}
	double P = func(X[0]);
	for (int i = 0;i < n - 1;i++) {
		P = P + coefs[i] * newt_calc(X, i, X_st);
	}
	cout << "P_3(X*) = " << P << '\n';;
	cout << "y(X*) = " << func(X_st) << '\n';
	cout << "Δ(P_3(X*)) = " << abs(P - func(X_st)) << '\n';
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	cout << "Interpolation for X_1 is: \n\n";
	interp_lagr(X_1, 4);
	interp_newt(X_1, 4);
	cout << "\n\n\n";
	cout << "Interpolation for X_2 is: \n\n";
	interp_lagr(X_2, 4);
	interp_newt(X_2, 4);
}
