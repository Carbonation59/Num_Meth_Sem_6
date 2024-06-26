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
	return 1 / (256 - x * x * x * x);
}

double X_0 = -2;
double X_k = 2;
double h_1 = 1;
double h_2 = 0.5;
double ans = 0.015827;

/*double func(double x) {
	return x / (3 * x + 4) / (3 * x + 4);
}

double X_0 = -1;
double X_k = 1;
double h_1 = 0.5;
double h_2 = 0.25;
double ans = -0.16474;*/

double rect_method(double h, int i) {
	double F = 0;
	for (double x = X_0; x < X_0 + i * h; x = x + h) {
		F = F + func((x + x + h) / 2);
	}
	F = F * h;
	return F;
}

double trapez_method(double h, int i) {
	double F = func(X_0) / 2;
	for (double x = X_0 + h; x <= X_0 + i * h; x = x + h) {
		F = F + func(x);
	}
	F = F - func(X_0 + h * i) / 2;
	F = F * h;
	return F;
}

double simpson_method(double h, int i) {
	double F = func(X_0);
	int k = 2;
	for (double x = X_0 + h; x <= X_0 + i * h; x = x + h) {
		k = (k >> 1);
		if (k == 1) {
			k = 4;
		}
		F = F + func(x) * k;
	}
	F = F - func(X_0 + h * i) * (k - 1);
	F = F * h / 3;
	return F;
}

double rrr_method(double F, double Fk, double p) {  // p - точность расчётной формулы
	return F + (F - Fk) / (pow(0.5, p) - 1);
}

void solve_intrgr() {
	cout.precision(7);
	cout << fixed;
	cout << "Integration with step = " << h_1 << ":\n\n";
	cout << "┌───────────┬───────────┬────────────┬────────────┬─────────────┬────────────┐\n";
	cout << "│     i     │    x_i    │    y_i     │ Rectangle  │ Trapezoidal │  Simpson   │\n";
	int k = 0;
	double x = X_0;
	while (x <= X_k) {
		cout << "│───────────┼───────────┼────────────┼────────────┼─────────────┼────────────│\n";
		cout << "│     " << k << "     │ ";
		if (x >= 0) {
			cout << '+';
		}
		cout << x << "  │  ";
		if (func(x) >= 0) {
			cout << '+';
		}
		cout << func(x) << "  │  ";
		if (rect_method(h_1, k) >= 0) {
			cout << '+';
		}
		cout << rect_method(h_1, k) << "  │  ";
		if (trapez_method(h_1, k) >= 0) {
			cout << '+';
		}
		cout << trapez_method(h_1, k) << "   │  ";
		if (simpson_method(h_1, k) >= 0) {
			cout << '+';
		}
		cout << simpson_method(h_1, k) << "  │\n";
		k++;
		x = x + h_1;
	}
	cout << "└───────────┴───────────┴────────────┴────────────┴─────────────┴────────────┘\n\n\n";
	cout << "Integration with step = " << h_2 << ":\n\n";
	cout << "┌───────────┬───────────┬────────────┬────────────┬─────────────┬────────────┐\n";
	cout << "│     i     │    x_i    │    y_i     │ Rectangle  │ Trapezoidal │  Simpson   │\n";
	k = 0;
	x = X_0;
	while (x <= X_k) {
		cout << "│───────────┼───────────┼────────────┼────────────┼─────────────┼────────────│\n";
		cout << "│     " << k << "     │ ";
		if (x >= 0) {
			cout << '+';
		}
		cout << x << "  │  ";
		if (func(x) >= 0) {
			cout << '+';
		}
		cout << func(x) << "  │  ";
		if (rect_method(h_2, k) >= 0) {
			cout << '+';
		}
		cout << rect_method(h_2, k) << "  │  ";
		if (trapez_method(h_2, k) >= 0) {
			cout << '+';
		}
		cout << trapez_method(h_2, k) << "   │  ";
		if (simpson_method(h_2, k) >= 0) {
			cout << '+';
		}
		cout << simpson_method(h_2, k) << "  │\n";
		k++;
		x = x + h_2;
	}
	cout << "└───────────┴───────────┴────────────┴────────────┴─────────────┴────────────┘\n\n";
	cout << "Exact value is: " << ans << "\n\n";
	cout << "Runge-Romberg-Richardson method:\n\n";
	k--;
	cout << "┌─────────────┬─────────────┬─────────────┐\n";
	cout << "│  Rectangle  │ Trapezoidal │   Simpson   │\n";
	cout << "│─────────────┼─────────────┼─────────────│\n";
	cout << "│  ";
	cout << abs(ans - rrr_method(rect_method(h_1, (k / 2 + k % 2)), rect_method(h_2, k), 2)) << "  │  ";
	cout << '+';
	cout << abs(ans - rrr_method(trapez_method(h_1, (k / 2 + k % 2)), trapez_method(h_2, k), 2)) << " │  ";
	cout << '+';
	cout << abs(ans - rrr_method(simpson_method(h_1, (k / 2 + k % 2)), simpson_method(h_2, k), 2)) << " │\n";
	cout << "│─────────────┼─────────────┼─────────────│\n";
	cout << "│  ";
	cout << '+';
	cout << abs(ans - rect_method(h_2, k)) << " │  ";
	cout << '+';
	cout << abs(ans - trapez_method(h_2, k)) << " │  ";
	cout << '+';
	cout << abs(ans - simpson_method(h_2, k)) << " │\n";
	cout << "└─────────────┴─────────────┴─────────────┘\n\n";
}


int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	solve_intrgr();
}
