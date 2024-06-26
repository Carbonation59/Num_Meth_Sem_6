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
	cout << "-6 * x1 + 5 * x2 = 51\n";
	cout << "-x1 + 13 * x2 + 6 * x3 = 100\n";
	cout << "-9 * x2 - 15 * x3 - 4 * x4 = -12\n";
	cout << "-x3 - 7 * x4 + x5 = 47\n";
	cout << "9 * x4 - 18 * x5 = -90\n";
	cout << '\n';
}

vector<double> sweep_method(vector<vector<double>> mtrx_A, vector<double> column_d) {
	int n = mtrx_A.size();
	vector<double> P(n);
	vector<double> Q(n);
	for (int i = 0;i < 5;i++) {
		if (i == 0) {
			P[i] = (-1) * mtrx_A[i][2] / mtrx_A[i][1];
			Q[i] = column_d[i] / mtrx_A[i][1];
		}
		else {
			P[i] = (-1) * mtrx_A[i][2] / (mtrx_A[i][1] + mtrx_A[i][0] * P[i - 1]);
			Q[i] = (column_d[i] - mtrx_A[i][0] * Q[i - 1]) / (mtrx_A[i][1] + mtrx_A[i][0] * P[i - 1]);
		}
	}
	vector<double> x(n);
	x[n - 1] = Q[n - 1];
	for (int i = n - 2;i > -1;i--) {
		x[i] = P[i] * x[i + 1] + Q[i];
	}
	return x;
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	vector<vector<double>> mtrx_A(5, vector<double>(3));
	print_task();
	mtrx_A[0][0] = 0; mtrx_A[0][1] = -6; mtrx_A[0][2] = 5;
	mtrx_A[1][0] = -1; mtrx_A[1][1] = 13; mtrx_A[1][2] = 6;
	mtrx_A[2][0] = -9; mtrx_A[2][1] = -15; mtrx_A[2][2] = -4;
	mtrx_A[3][0] = -1; mtrx_A[3][1] = -7; mtrx_A[3][2] = 1;
	mtrx_A[4][0] = 9; mtrx_A[4][1] = -18; mtrx_A[4][2] = 0;
	vector<double> column_d(5);
	column_d[0] = 51; column_d[1] = 100; column_d[2] = -12; column_d[3] = 47; column_d[4] = -90;
	vector<double> x = sweep_method(mtrx_A, column_d);
	cout << "Column x is:\n\n";
	for (int i = 0;i < 5;i++) {
		cout << x[i] << '\n';
	}

}
