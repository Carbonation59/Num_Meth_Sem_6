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

double X_0 = 1;
double X_k = 2;
double h1 = 0.1;
double Y_0 = exp(1) + 2;
double Z_0 = exp(1) + 1;

double func(double x) {
	return x + 1 + exp(x);
}

double func_f(double x, double y, double z) {
	return z;
}

double func_g(double x, double y, double z) {
	return ((x + 1) * z - y) / x;
}

void eiler_method() {
	cout << "Eiler method:\n\n";
	cout.precision(8);
	cout << fixed;
	double h = h1;
	double cur_y = Y_0;
	double cur_z = Z_0;
	double delta_y;
	double delta_z;
	double app;
	int k = 0;
	cout << "┌───────────┬───────────┬────────────┬────────────┬─────────────┬────────────┬────────────┬──────────────┬─────────────┐\n\n";
	cout << "│     k     │    x_k    │    y_k     │     z_k    │    Δy_k     │    Δz_k    │   check    │    y_true    │    eps_k    │\n\n";
	for (double x = X_0; x < X_k; x = x + h) {
		cout << "│───────────┼───────────┼────────────┼────────────┼─────────────┼────────────┼────────────┼──────────────┼─────────────│\n\n";
		cout << "│     " << k;
		if (k < 10) {
			cout << ' ';
		}
		app = abs((pow(cur_y, h) - pow(cur_y, 2 * h)) / 3);  // 2^2 - 1 = 3
		delta_y = h * func_f(x, cur_y, cur_z);
		delta_z = h * func_g(x, cur_y, cur_z);
		cout << "    │ " << x << "│ " << cur_y;
		if (cur_y < 10) {
			cout << ' ';
		}
		cout << "│ " << cur_z << " │  " << delta_y << " │ " << delta_z << " │ " << app << " │  " << func(x) << "  │ " << abs(func(x) - cur_y) << "  │\n\n";
		if (app > abs(func(x) - cur_y)) {
			h = h / 2;
		}
		cur_y = cur_y + delta_y;
		cur_z = cur_z + delta_z;
		k++;
	}
	cout << "│───────────┼───────────┼────────────┼────────────┼─────────────┼────────────┼────────────┼──────────────┼─────────────│\n\n";
	cout << "│    " << k << "     │ ";
	cout << X_k << "│ " << cur_y << "│ " << cur_z << " │             │            │            │  " << func(X_k) << " │ " << abs(func(X_k) - cur_y) << "  │\n\n";
	cout << "└───────────┴───────────┴────────────┴────────────┴─────────────┴────────────┴────────────┴──────────────┴─────────────┘\n\n";
	cout << "Columns x_k ans y_k - solve\n\n\n";
}

void runge_kutt_method() {
	cout << "Runge-Kutt method:\n\n";
	cout.precision(8);
	cout << fixed;
	double h = h1;
	double cur_y = Y_0;
	double cur_z = Z_0;
	double delta_y;
	double delta_z;
	double app;
	int k = 0;
	cout << "┌───────────┬───────────┬────────────┬────────────┬─────────────┬────────────┬────────────┬──────────────┬─────────────┐\n\n";
	cout << "│     k     │    x_k    │    y_k     │     z_k    │    Δy_k     │    Δz_k    │   check    │    y_true    │    eps_k    │\n\n";
	for (double x = X_0; x < X_k; x = x + h) {
		cout << "│───────────┼───────────┼────────────┼────────────┼─────────────┼────────────┼────────────┼──────────────┼─────────────│\n\n";
		vector<double> K(4);
		vector<double> L(4);
		K[0] = h * func_f(x, cur_y, cur_z);
		L[0] = h * func_g(x, cur_y, cur_z);
		for (int j = 1;j < 3;j++) {
			K[j] = h * func_f(x + h / 2, cur_y + K[j - 1] / 2, cur_z + L[j - 1] / 2);
			L[j] = h * func_g(x + h / 2, cur_y + K[j - 1] / 2, cur_z + L[j - 1] / 2);
		}
		K[3] = h * func_f(x + h, cur_y + K[2], cur_z + L[2]);
		L[3] = h * func_g(x + h, cur_y + K[2], cur_z + L[2]);
		app = max(abs((K[1] - K[2]) / (K[0] - K[1])), abs((L[1] - L[2]) / (L[0] - L[1])));
		delta_y = (K[0] + 2 * K[1] + 2 * K[2] + K[3]) / 6;
		delta_z = (L[0] + 2 * L[1] + 2 * L[2] + L[3]) / 6;
		cout << "│     " << k << "     │ ";
		cout << x << "│ " << cur_y << " │ " << cur_z << " │  " << delta_y << " │ " << delta_z << " │ " << app << " │  " << func(x) << "  │ " << abs(func(x) - cur_y) << "  │\n\n";
		if (app > 0.1) {
			h = h / 2;
		}
		if (app < 0.01) {
			h = h * 2;
		}
		cur_y = cur_y + delta_y;
		cur_z = cur_z + delta_z;
		k++;
	}
	cout << "│───────────┼───────────┼────────────┼────────────┼─────────────┼────────────┼────────────┼──────────────┼─────────────│\n\n";
	cout << "│    " << k << "     │ ";
	cout << X_k << "│ " << cur_y << "│ " << cur_z << " │             │            │            │  " << func(X_k) << " │ " << abs(func(X_k) - cur_y) << "  │\n\n";
	cout << "└───────────┴───────────┴────────────┴────────────┴─────────────┴────────────┴────────────┴──────────────┴─────────────┘\n\n";
	cout << "Columns x_k ans y_k - solve\n\n\n";
}

void adams_method() {
	cout << "Adams method:\n\n";
	cout.precision(8);
	cout << fixed;
	double h = h1;
	double cur_y = Y_0;
	double cur_z = Z_0;
	double delta_y;
	double delta_z;
	vector<double> arr_y;
	vector<double> arr_f;
	vector<double> arr_z;
	vector<double> arr_g;
	int k = 0;
	cout << "┌───────────┬───────────┬────────────┬────────────┬─────────────────────────┬──────────────┬─────────────┐\n";
	cout << "│     k     │    x_k    │    y_k     │     z_k    │     f(x_k, y_k, z_k)    │    y_true    │    eps_k    │\n";
	for (double x = X_0; x < X_k; x = x + h) {
		cout << "│───────────┼───────────┼────────────┼────────────┼─────────────────────────┼──────────────┼─────────────│\n";
		cout << "│     " << k << "     │ ";
		cout << x << "│ " << cur_y << " │ " << cur_z << " │        " << func_f(x, cur_y, cur_z) << "       │  " << func(x) << "  │ " << abs(func(x) - cur_y) << "  │\n";
		if (k < 4) {
			vector<double> K(4);
			vector<double> L(4);
			K[0] = h * func_f(x, cur_y, cur_z);
			L[0] = h * func_g(x, cur_y, cur_z);
			for (int j = 1;j < 3;j++) {
				K[j] = h * func_f(x + h / 2, cur_y + K[j - 1] / 2, cur_z + L[j - 1] / 2);
				L[j] = h * func_g(x + h / 2, cur_y + K[j - 1] / 2, cur_z + L[j - 1] / 2);
			}
			K[3] = h * func_f(x + h, cur_y + K[2], cur_z + L[2]);
			L[3] = h * func_g(x + h, cur_y + K[2], cur_z + L[2]);
			delta_z = (L[0] + 2 * L[1] + 2 * L[2] + L[3]) / 6;
			cur_z = cur_z + delta_z;
			delta_y = (K[0] + 2 * K[1] + 2 * K[2] + K[3]) / 6;
			cur_y = cur_y + delta_y;
		}
		else {
			cur_y = arr_y[k - 1] + h * (55 * arr_f[k - 1] - 59 * arr_f[k - 2] + 37 * arr_f[k - 3] - 9 * arr_f[k - 4]) / 24;
			cur_z = arr_z[k - 1] + h * (55 * arr_g[k - 1] - 59 * arr_g[k - 2] + 37 * arr_g[k - 3] - 9 * arr_g[k - 4]) / 24;
		}
		arr_y.push_back(cur_y);
		arr_f.push_back(func_f(x, cur_y, cur_z));
		arr_z.push_back(cur_z);
		arr_g.push_back(func_g(x, cur_y, cur_z));
		k++;
	}
	cout << "│───────────┼───────────┼────────────┼────────────┼─────────────────────────┼──────────────┼─────────────│\n";
	cout << "│    " << k << "     │ ";
	cout << X_k << "│ " << cur_y << "│ " << cur_z << " │        " << func_f(X_k, cur_y, cur_z) << "       │  " << func(X_k) << " │ " << abs(func(X_k) - cur_y) << "  │\n";
	cout << "└───────────┴───────────┴────────────┴────────────┴─────────────────────────┴──────────────┴─────────────┘\n\n";
	cout << "Columns x_k ans y_k - solve\n\n\n";
}


int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	eiler_method();
	runge_kutt_method();
	adams_method();
}
