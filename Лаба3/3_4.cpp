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

double X[5] = { -0.2, 0, 0.2, 0.4, 0.6 };
double Y[5] = { -0.40136, 0, 0.40136, 0.81152, 1.2435 };

double X_st = 0.2;

//double X[5] = { 0, 0.1, 0.2, 0.3, 0.4 };
//double Y[5] = { 1, 1.1052, 1.2214, 1.3499, 1.4918 };

//double X_st = 0.2;

pair<double, double> calc(int j) {
	double y1 = (Y[j] - Y[j - 1]) / (X[j] - X[j - 1]);
	double y2 = (Y[j + 1] - Y[j]) / (X[j + 1] - X[j]);
	cout.precision(5);
	cout << fixed;
	cout << "Left-Hand:" << y1 << '\n';
	cout << "Right-Hand:" << y2 << "\n";
	cout << "Half_Sum: " << (y1 + y2) / 2 << '\n';  
	double y3 = y1 + (y2 - y1) / (X[j + 1] - X[j - 1]) * (2 * X_st - X[j - 1] - X[j]);
	cout << "Second order accuracy: " << y3 << "\n\n\n"; // если совпадают, то сетка равномерна
	return { y1, y2 };
}

void fst_scd_der(int n) {
	cout << "First derivative is:\n\n";
	pair<double, double> y;
	int j;
	for (int i = 0;i < n;i++) {
		if (X_st == X[i]) {
			y = calc(i);
			j = i;
		}
	}
	cout << "Second derivative is: ";
	cout << 2 * (y.second - y.first) / (X[j + 1] - X[j - 1]) << '\n';
}


int main() {
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	fst_scd_der(5);
}
