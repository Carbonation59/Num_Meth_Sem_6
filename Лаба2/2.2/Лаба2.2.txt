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

double func1(double x1, double x2) {
	return 3 * x1 - cos(x2);
}

double d_func1_d_x1(double x1, double x2) {
	return 3;
}

double d_func1_d_x2(double x1, double x2) {
	return sin(x2);
}

double func2(double x1, double x2) {
	return 3 * x2 - exp(x1);
}

double d_func2_d_x1(double x1, double x2) {
	return - exp(x1);
}

double d_func2_d_x2(double x1, double x2) {
	return 3;
}

double det_A1(double x1, double x2) {
	return func1(x1, x2) * d_func2_d_x2(x1, x2) - func2(x1, x2) * d_func1_d_x2(x1, x2);
}

double det_A2(double x1, double x2) {
	return func2(x1, x2) * d_func1_d_x1(x1, x2) - func1(x1, x2) * d_func2_d_x1(x1, x2);
}

double det_J(double x1, double x2) {
	return d_func1_d_x1(x1, x2) * d_func2_d_x2(x1, x2) - d_func1_d_x2(x1, x2) * d_func2_d_x1(x1, x2);
}


void newton_method(double a1, double b1, double a2, double b2, double eps) {
	double x11 = (a1 + b1) / 2;
	double x21 = (a2 + b2) / 2;
	cout << "Newton method:\n\n";
	cout << "┌─────────┬────────────┬────────────────────┬─────────────────────────────┬─────────────────────────────┬────────────┬────────────┬────────────┐\n";
	cout << "│    k    │    x1_k    │   f1(x1_k, x2_k)   │   d_f1(x1_k, x2_k) / d_x1   │   d_f1(x1_k, x2_k) / d_x2   │ det_A1(k)  │ det_A2(k)  │  det_J(k)  │\n";
	cout << "│         │    x2_k    │   f2(x1_k, x2_k)   │   d_f2(x1_k, x2_k) / d_x1   │   d_f2(x1_k, x2_k) / d_x2   │            │            │            │\n";
	cout << "│─────────┼────────────┼────────────────────┼─────────────────────────────┼─────────────────────────────┼────────────┼────────────┼────────────┼\n";
	cout.precision(5);
	cout << fixed;
	double x12 = x11 - det_A1(x11, x21) / det_J(x11, x21);
	double x22 = x21 - det_A2(x11, x21) / det_J(x11, x21);
	int k = 0;
	cout << "│    " << k << "    " << "│   " << x11 << "  " << "│      " << func1(x11, x21) << "      " << "│            " << d_func1_d_x1(x11, x21) << "          " << "│           " << d_func1_d_x2(x11, x21)<< "           " << "│  " << det_A1(x11, x21) << "  " << "│  "<< det_A2(x11, x21) << "  " << "│   " << det_J(x11, x21) << "  │\n";
	cout << "│         " << "│   " << x21 << "  " << "│      " << func2(x11, x21) << "      " << "│           " << d_func2_d_x1(x11, x21) << "          " << "│           " << d_func2_d_x2(x11, x21) << "           " << "│            " << "│            " << "│            │\n";
	cout << "│─────────┼────────────┼────────────────────┼─────────────────────────────┼─────────────────────────────┼────────────┼────────────┼────────────┼\n";
	k++;
	double mx = max(abs(x22 - x21), abs(x12 - x11));
	while (mx >= eps) {
		x11 = x12;
		x21 = x22;
		cout << "│    " << k << "    " << "│   " << x11 << "  " << "│      " << func1(x11, x21) << "      " << "│            " << d_func1_d_x1(x11, x21) << "          " << "│           " << d_func1_d_x2(x11, x21) << "           " << "│  " << det_A1(x11, x21) << "  " << "│  " << det_A2(x11, x21) << "  " << "│   " << det_J(x11, x21) << "  │\n";
		cout << "│         " << "│   " << x21 << "  " << "│      " << func2(x11, x21) << "      " << "│           " << d_func2_d_x1(x11, x21) << "          " << "│           " << d_func2_d_x2(x11, x21) << "           " << "│            " << "│            " << "│            │\n";
		cout << "│─────────┼────────────┼────────────────────┼─────────────────────────────┼─────────────────────────────┼────────────┼────────────┼────────────┼\n";
		x12 = x11 - det_A1(x11, x21) / det_J(x11, x21);
		x22 = x21 - det_A2(x11, x21) / det_J(x11, x21);
		k++;
		mx = max(abs(x22 - x21), abs(x12 - x11));
	}
	x11 = x12;
	x21 = x22;
	cout << "│    " << k << "    " << "│   " << x11 << "  " << "│\n";
	cout << "│         " << "│   " << x21 << "  │\n";
	cout << "└─────────┴────────────┘\n\n";
	cout << "x1* ≈ " << x11 << ", x2* ≈ " << x21 << "\n\n";
}

double phi1(double x1, double x2) {
	return cos(x2) / 3;
}

double d_phi1_d_x1(double x1, double x2) {
	return 0;
}

double d_phi1_d_x2(double x1, double x2) {
	return -sin(x2) / 3;
}

double phi2(double x1, double x2) {
	return exp(x1) / 3;
}

double d_phi2_d_x1(double x1, double x2) {
	return exp(x1) / 3;
}

double d_phi2_d_x2(double x1, double x2) {
	return 0;
}

double find_q(double a1, double b1, double a2, double b2) {
	double mx = 0;
	for (double i = a1; i <= b1;i = i + 0.001) {
		for (double j = a2; j <= b2;j = j + 0.001) {
			mx = max(mx, abs(d_phi1_d_x1(i, j)) + abs(d_phi1_d_x2(i, j)));
		}
	}
	for (double i = a1; i <= b1;i = i + 0.001) {
		for (double j = a2; j <= b2;j = j + 0.001) {
			mx = max(mx, abs(d_phi2_d_x1(i, j)) + abs(d_phi2_d_x2(i, j)));
		}
	}
	return mx;
}

void method_simp_iter(double a1, double b1, double a2, double b2, double eps) {
	double q = find_q(a1, b1, a2, b2);
	double x11 = (a1 + b1) / 2;
	double x21 = (a2 + b2) / 2;
	cout << "Method of simple iterations:\n\n";
	if (q >= 1) {
		cout << "Error: wrong fucntion phi\n\n";
		return;
	}
	cout << "q is: " << q << "\n\n";
	cout << "┌───────────────┬─────────────┬──────────────────────┐\n";
	cout << "│       k       │    x1_k     │   phi1(x1_k, x2_k)   │\n";
	cout << "│               │    x2_k     │   phi2(x1_k, x2_k)   │\n";
	cout << "│───────────────┼─────────────┼──────────────────────|\n";
	cout.precision(5);
	cout << fixed;
	double x12 = phi1(x11, x21);
	double x22 = phi2(x11, x21);
	int k = 0;
	cout << "│      " << k << "        " << "│   " << x11 << "   " << "│       " << x12 << "        │\n";
	cout << "│               " << "│   " << x21 << "   " << "│       " << x22 << "        │\n";
	cout << "│───────────────┼─────────────┼──────────────────────|\n";
	k++;
	double mx = max(abs(x22 - x21), abs(x12 - x11));
	while (q / (1 - q) * mx >= eps) {
		x11 = x12;
		x21 = x22;
		x12 = phi1(x11, x21);
		x22 = phi2(x11, x21);
		cout << "│      " << k << "        " << "│   " << x11 << "   " << "│       " << x12 << "        │\n";
		cout << "│               " << "│   " << x21 << "   " << "│       " << x22 << "        │\n";
		cout << "│───────────────┼─────────────┼──────────────────────|\n";
		k++;
		mx = max(abs(x22 - x21), abs(x12 - x11));
	}
	x11 = x12;
	x21 = x22;
	cout << "│      " << k << "        " << "│   " << x11 << "   │\n";
	cout << "│               " << "│   " << x21 << "   │\n";
	cout << "└───────────────┴─────────────┘\n\n";
	cout << "x1* ≈ " << x11 << ", x2* ≈ " << x21 << "\n\n";

}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	double eps = 0.0001;
	newton_method(0, 0.5, 0, 0.5, eps);
	method_simp_iter(0, 0.5, 0, 0.5, eps);
}
