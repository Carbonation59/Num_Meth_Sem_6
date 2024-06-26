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

double func(double x) {
	return x - 3 + 1 / (x + 1);
}

double func_p(double x) {
	return x - 3;
}

double func_q(double x) {
	return -1;
}

double func_f(double x) {
	return 0;
}

double rrr_method(double F, double Fk, double p) {  // p - точность расчётной формулы
	return F + (F - Fk) / (pow(0.5, p) - 1);
}


double h = 0.2;
double X_0 = 0;
double X_n = 1;

vector<double> find_c(vector<vector<double>> matrix, vector<double> ans) {
	int n = matrix.size();
	double coef;
	for (int i = 0;i < n - 1;i++) {
		coef = matrix[i + 1][0] / matrix[i][1];
		for (int j = 1;j < 3;j++) {
			matrix[i][j] = matrix[i][j] * coef;
			matrix[i + 1][j - 1] = matrix[i + 1][j - 1] - matrix[i][j];
		}
		ans[i] = ans[i] * coef;
		ans[i + 1] = ans[i + 1] - ans[i];
	}
	vector<double> coefs(n);
	for (int i = n - 1;i >= 0;i--) {
		coefs[i] = ans[i] / matrix[i][1];
		if (i != 0) {
			ans[i - 1] = ans[i - 1] - coefs[i] * matrix[i - 1][2];
		}
	}
	return coefs;
}

void finite_difference_method(int n) {
	vector<vector<double>> matrix(n, vector<double>(3));
	vector<double> ans(n);
	double y_a = func(X_0);
	double y_b = func(X_n);
	matrix[0][1] = (-2 + h * h * func_q(X_0 + h));
	matrix[0][2] = (1 + h * func_p(X_0 + h) / 2);
	ans[0] = h * h * func_f(X_0 + h) - (1 - h * func_p(X_0 + h) / 2) * y_a;
	for (int i = 1;i < n - 1;i++) {
		matrix[i][0] = (1 - h * func_p(X_0 + h + h * i) / 2);
		matrix[i][1] = (-2 + h * h * func_q(X_0 + h + h * i));
		matrix[i][2] = (1 + h * func_p(X_0 + h + h * i) / 2);
		ans[i] = h * h * func_f(X_0 + h + h * i);
	}
	matrix[n - 1][0] = (1 - h * func_p(X_0 + h * n) / 2);
	matrix[n - 1][1] = (-2 + h * h * func_q(X_0 + h * n));
	ans[n - 1] = h * h * func_f(X_0 + h * n) - (1 + h * func_p(X_0 + h * n) / 2) * y_b;
	vector<double> y_n = find_c(matrix, ans);
	cout.precision(5);
	cout << fixed;
	cout << "┌───────────┬──────────┬───────────┬───────────┬───────────┬───────────┬───────────┐\n";
	cout << "│     k     │     0    │     1     │     2     │     3     │     4     │     5     │\n";
	cout << "│───────────┼──────────┼───────────┼───────────┼───────────┼───────────┼───────────│\n";
	cout << "│  " << "  x_k    │ " << X_0 << "  │  " << X_0 + h << "  │  " << X_0 + h * 2 << "  │  " << X_0 + h * 3 << "  │  " << X_0 + h * 4 << "  │  " << X_0 + h * 5 << "  │\n";
	cout << "│───────────┼──────────┼───────────┼───────────┼───────────┼───────────┼───────────│\n";
	cout << "│  " << "  y_k    │ " << y_a << " │  " << y_n[0] << " │  " << y_n[1] << " │  " << y_n[2] << " │  " << y_n[3] << " │ " << y_n[4] << "  │\n";
	cout << "│───────────┼──────────┼───────────┼───────────┼───────────┼───────────┼───────────│\n";
	cout << "│  " << "Run_Rom  │ " << func(X_0) - rrr_method(y_a, func(X_0), 4) << "  │  " << func(X_0 + h) - rrr_method(y_n[0], func(X_0 + h), 4) << "  │  " << func(X_0 + h * 2) - rrr_method(y_n[1], func(X_0 + h * 2), 4) << "  │  " << func(X_0 + h * 3) - rrr_method(y_n[2], func(X_0 + h * 3), 4) << "  │  " << func(X_0 + h * 4) - rrr_method(y_n[3], func(X_0 + h * 4), 4) << "  │ " << func(X_0 + h * 5) - rrr_method(y_n[4], func(X_0 + h * 5), 4) << "  │\n";
	cout << "└───────────┴──────────┴───────────┴───────────┴───────────┴───────────┴───────────┘\n\n";
}

//Shooting

double Func(double x, double y, double z) {
	return exp(x) + sin(y);
	//	return ((x+1)*z+2*(x-1)*y)/x;
}

double Yist(double x) {
	//	return exp(x)+sin(y);
	return exp(2.0 * x) + (3.0 * x + 1.0) * exp(-x);
}

double eps = 0.001;

vector<double> RRR(vector<double> x, vector<double> y, double p) {
	vector<double> z;
	for (int i = 0; i < x.size() - 1; i++) {
		z.push_back(abs((x[i] - y[2 * i]) / (pow(2, p) - 1)));
	}
	return z;
}

pair<vector<double>, vector<double>> RK(double h, double a, double b, double y0, double z0) {
	vector<double> x, y, z;
	double i, k1, k2, k3, l1, l2, l3, k4, l4;
	x.push_back(a);
	y.push_back(y0);
	z.push_back(z0);
	int k = 0;
	for (i = a + h; i <= b + h; i = i + h) {
		k1 = h * z[k];
		l1 = h * Func(x[k], y[k], z[k]);
		k2 = h * (z[k] + 1.0 / 2 * l1);
		l2 = h * Func(x[k] + 1.0 / 2 * h, y[k] + 1.0 / 2 * k1, z[k] + 1.0 / 2 * l1);
		k3 = h * (z[k] + 1.0 / 2 * l2);
		l3 = h * Func(x[k] + 1.0 / 2 * h, y[k] + 1.0 / 2 * k2, z[k] + 1.0 / 2 * l2);
		k4 = h * (z[k] + l3);
		l4 = h * Func(x[k] + h, y[k] + k3, z[k] + l3);

		x.push_back(i);
		y.push_back(y[k] + 1.0 / 6 * (k1 + 2 * k2 + 2 * k3 + k4));
		z.push_back(z[k] + 1.0 / 6 * (l1 + 2 * l2 + 2 * l3 + l4));
		k++;
	}
	return { y,z };
}


vector<double> shoot(double h, double a, double b) {
	double et_prev = 1.0;
	double et_i = 0.8;
	pair<vector<double>, vector<double>> tmp;
	tmp = RK(h, a, b, 1.0, et_prev);
	vector<double> y_prev = tmp.first;
	vector<double> z_prev = tmp.second;
	tmp = RK(h, a, b, 1.0, et_i);
	vector<double> y_i = tmp.first;
	vector<double> z_i = tmp.second;
	double Fi_prev = y_prev[y_prev.size() - 2] - 2.0;
	double Fi_i = y_i[y_i.size() - 2] - 2.0;
	double t;
	while (abs(Fi_i) > eps) {
		t = et_i;
		et_i = et_i - Fi_i * (et_i - et_prev) / (Fi_i - Fi_prev);
		et_prev = t;
		y_prev = y_i;
		tmp = RK(h, a, b, 1, et_i);
		y_i = tmp.first;
		z_i = tmp.second;
		Fi_prev = Fi_i;
		Fi_i = y_i[y_i.size() - 2] - 2.0;

	}
	return y_i;
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	finite_difference_method(5);
	double h = 0.1;
	double a, b;
	a = 0;
	b = 1;
	vector<double> x, x1, x2;

	x = shoot(h, a, b);
	x1 = shoot(h / 2.0, a, b);
	x2 = RRR(x, x1, 4);
	//cout << "Метод Рунге-Кутты 4-ого порядка:\n";
	cout << "Shooting method:\n\n";
	cout << "x y y(ист) ek e(RRR)\n";
	for (int i = 0; i < x.size() - 1; i++) {
		cout << a + i * h << " " << x[i] << " " << Yist(a + i * h) << " " << abs(Yist(a + i * h) - x[i]) << " " << x2[i] << '\n';//	
	}
	cout << '\n';
}
