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
	return pow(4, x) - 5 * x - 2;
}

double func_1(double x) {
	return log(4) * pow(4, x) - 5;
}


double func_2(double x) {
	return log(4) * log(4) * pow(4, x);
}

void newton_method(double x, double eps) {
	cout << "Newton method:\n\n";
	cout << "f(x) * f''(x) = " << func(x) * func_2(x) << " > 0 \n\n";
	cout << "┌───────────────┬─────────────┬────────────┬────────────┬──────────────────────┐\n";
	cout << "│      k        │     x_k     │   f(x_k)   │   f'(x_k)  │  - f(x_k) / f'(x_k)  │\n";
	cout << "│───────────────┼─────────────┼────────────┼────────────┼──────────────────────|\n";
	cout.precision(4);
	cout << fixed;
	double x1 = x - func(x) / func_1(x);
	int k = 0;
	cout << "│      " << k << "        ";
	cout << "│   " << x << "    ";
	cout << "│   " << func(x) << "   ";
	cout << "|";
	if (func_1(x) < 10) {
		cout << " ";
	}
	cout << "   " << func_1(x) << "  ";
	cout << "│        " << -func(x) / func_1(x) << "       │\n";
	k++;
	while (abs(x1 - x) >= eps) {
		cout << "│      " << k << "        ";
		cout << "│   " << x1 << "    ";
		cout << "│   " << func(x1) << "   ";
		cout << "|";
		if (func_1(x1) < 10) {
			cout << " ";
		}
		cout << "   " << func_1(x1) << "  ";
		cout << "│        " << -func(x1) / func_1(x1) << "       │\n";
		x = x1;
		x1 = x - func(x) / func_1(x);
		k++;
	}
	cout << "│      " << k << "        ";
	cout << "│  [" << x1 << "]   ";
	cout << "│            ";
	cout << "│            ";
	cout << "│                      │\n";
	cout << "└───────────────┴─────────────┴────────────┴────────────┴──────────────────────┘\n\n";
	cout << "x* ≈ " << x1 << "\n\n";
}

double get_sign(double x) {
	if (x < 0) {
		return -1;
	}
	return 1;
}

double mx_func_1 = 11.81;

double func_phi(double x) {     // монотонна, поэтому решение точно найдётся
	return x - func(x) * get_sign(func_1(x)) / mx_func_1;
}

double func_phi_1(double x) {     // монотонна, поэтому решение точно найдётся
	return 1 - func_1(x) * get_sign(func_1(x)) / mx_func_1;
}

bool check_phi(double a, double b) {
	double mx = -100;
	double mn = 100;
	for (double i = a; i <= b; i = i + 0.01) {
		mx = max(mx, func_phi(i));
		mn = min(mn, func_phi(i));
	}
	if (mx <= b && mn >= a) {
		return true;
	}
	return false;
}

double find_q(double a, double b) {
	double mx = 0;
	for (double i = a; i <= b; i = i + 0.01) {
		mx = max(mx, abs(func_phi_1(i)));
	}
	return mx;
}

void method_simp_iter(double a, double b, double eps) {
	double q = find_q(a, b);
	double x = (a + b) / 2;
	bool test = check_phi(a, b);
	cout << "Method of simple iterations :\n\n";
	if (!test) {
		cout << "Error^ wrong fucntion phi\n\n";
		return;
	}
	cout << "phi(x) in [a, b]: true\n";
	cout << "q is: " << q << "\n\n";
	cout << "┌───────────────┬─────────────┬──────────────┐\n";
	cout << "│      k        │     x_k     │   phi(x_k)   │\n";
	cout << "│───────────────┼─────────────┼──────────────|\n";
	cout.precision(4);
	cout << fixed;
	double x1 = func_phi(x);
	int k = 0;
	cout << "│      " << k << "        ";
	cout << "│   " << x << "    ";
	cout << "│    " << func_phi(x) << "    |\n";
	k++;
	while (q / (1 - q) * abs(x1 - x) > eps) {
		cout << "│      " << k << "        ";
		cout << "│   " << x1 << "    ";
		cout << "│    " << func_phi(x1) << "    |\n";
		x = x1;
		x1 = func_phi(x);
		k++;
	}
	cout << "│      " << k << "        ";
	cout << "│  [" << x1 << "]   ";
	cout << "│              │\n";
	cout << "└───────────────┴─────────────┴──────────────┘\n\n";
	cout << "x* ≈ " << x1 << '\n';

}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	double eps = 0.001;
	newton_method(1.8, eps);
	method_simp_iter(1.6, 1.8, eps);
}
