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

double X_st = -0.5;

double X[5] = { -2, -1, 0, 1, 2 };
double func[5] = { -1.8647, -0.63212, 1, 3.7183, 9.3891 };

vector<double> find_c(vector<double> h, int n) {
	vector<vector<double>> matrix;
	vector<double> ans;
	double c1, c2, c3, cur_ans;
	c2 = 2 * (h[1] + h[2]);
	c3 = h[2];
	matrix.push_back({ 0, c2, c3 });
	cur_ans = 3 * ((func[2] - func[1]) / h[2] - (func[1] - func[0]) / h[1]);
	ans.push_back(cur_ans);
	for (int i = 3;i <= n - 1;i++) {
		c1 = h[i - 1];
		c2 = 2 * (h[i - 1] + h[i]);
		c3 = h[i];
		matrix.push_back({ c1, c2, c3 });
		cur_ans = 3 * ((func[i] - func[i - 1]) / h[i] - (func[i - 1] - func[i - 2]) / h[i - 1]);
		ans.push_back(cur_ans);
	}
	c1 = h[n - 1];
	c2 = 2 * (h[n - 1] + h[n]);
	matrix.push_back({ c1, c2, 0 });
	cur_ans = 3 * ((func[n] - func[n - 1]) / h[n] - (func[n - 1] - func[n - 2]) / h[n - 1]);
	ans.push_back(cur_ans);
	double coef;
	for (int i = 0;i < n - 2;i++) {
		coef = matrix[i + 1][0] / matrix[i][1];
		for (int j = 1;j < n - 1;j++) {
			matrix[i][j] = matrix[i][j] * coef;
			matrix[i + 1][j - 1] = matrix[i + 1][j - 1] - matrix[i][j];
		}
		ans[i] = ans[i] * coef;
		ans[i + 1] = ans[i + 1] - ans[i];
	}
	vector<double> coefs(4);
	for (int i = n - 1;i > 0;i--) {
		coefs[i] = ans[i - 1] / matrix[i - 1][1];
		if (i != 1) {
			ans[i - 2] = ans[i - 2] - coefs[i] * matrix[i - 2][2];
		}
	}
	coefs[0] = 0;
	return coefs;
}

double calc(double x, double xi, double a, double b, double c, double d) {
	return a + (x - xi) * b + c * (x - xi) * (x - xi) + d * (x - xi) * (x - xi) * (x - xi);
}

void cub_sp(int n) {
	cout << "Cubic spline:\n\n";
	vector<double> h;
	h.push_back(0);
	for (int i = 1;i <= n;i++) {
		h.push_back(X[i] - X[i - 1]);
	}
	vector<double> coef_c = find_c(h, n);
	vector<double> coef_a(n);
	for (int i = 0;i < n;i++) {
		coef_a[i] = func[i];
	}
	vector<double> coef_b(n);
	for (int i = 0;i < n - 1;i++) {
		coef_b[i] = (func[i + 1] - func[i]) / h[i + 1] - h[i + 1] * (coef_c[i + 1] + 2 * coef_c[i]) / 3;
	}
	coef_b[n - 1] = (func[n] - func[n - 1]) / h[n] - 2 * h[n] * coef_c[n - 1] / 3;
	vector<double> coef_d(n);
	for (int i = 0;i < n - 1;i++) {
		coef_d[i] = (coef_c[i + 1] - coef_c[i]) / 3 / h[i + 1];
	}
	coef_d[n - 1] = (coef_c[n - 1] / 3 / h[n]) * (-1);
	cout.precision(5);
	cout << fixed;
	cout << "┌───────────────┬────────────────────┬────────────┬────────────┬────────────┬────────────┐\n";
	cout << "│       i       │   [x_(i-1), x_i]   │     a_i    │     b_i    │     c_i    │     d_i    │\n";
	for (int i = 0;i < n;i++) {
		cout << "│───────────────┼────────────────────┼────────────┼────────────┼────────────┼────────────│\n";
		cout << "│       " << i + 1 << "       │[";
		if (X[i] >= 0) {
			cout << '+';
		}
		cout << X[i] << "; ";
		if (X[i + 1] >= 0) {
			cout << '+';
		}
		cout << X[i + 1] << "]│  ";
		if (coef_a[i] >= 0) {
			cout << '+';
		}
		cout << coef_a[i] << "  │  ";
		if (coef_b[i] >= 0) {
			cout << '+';
		}
		cout << coef_b[i] << "  │  ";
		if (coef_c[i] >= 0) {
			cout << '+';
		}
		cout << coef_c[i] << "  │  ";
		if (coef_d[i] >= 0) {
			cout << '+';
		}
		cout << coef_d[i] << "  │   \n";
	}
	cout << "└───────────────┴────────────────────┴────────────┴────────────┴────────────┴────────────┘\n\n";
	for (int i = 0;i < n - 1;i++) {
		if (X[i] <= X_st && X_st <= X[i + 1]) {
			cout << "f(X*) = " << calc(X_st, X[i], coef_a[i], coef_b[i], coef_c[i], coef_d[i]) << '\n';
		}
	}
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	cub_sp(4);
}
