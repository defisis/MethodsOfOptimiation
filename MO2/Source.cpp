#include <iostream>
#include <vector>
#include <fstream>
#include "Vector.h"
#include <array>
#include "matrix.h"
#include <string>
#include <iomanip>

const double EPS = 1e-20;
int fcount = 0;

double func(Vector<double> x);
Vector<double> grad(Vector <double>);
double funcForSearch(double l, Vector<double> xc, matrix etta, Vector<double> gr);
void searchMinimumInterval(double & begin, double & end, double x0, double delta, Vector<double> xc, matrix etta, Vector<double> gr);
double goldenSection(double begin, double end, double eps, Vector<double> xc, matrix etta, Vector<double> gr);
void broidenMethod(Vector<double> &x, double eps, int maxiter);

int main() {
	Vector<double> x(2);
	std::cin >> x;
	double eps = 1e-12;
	int maxiter = 10000;
	broidenMethod(x, eps, maxiter);
}

double func(Vector<double> x) {
	return -2 * exp(-pow(x.items[0] - 3, 2) - pow(x.items[1] - 2, 2)) - 3 * exp(-pow(x.items[0] - 1, 2) - pow(x.items[1] - 1, 2) / 9);

	//return 4 * x.items[0] * x.items[0] + x.items[1] * x.items[1] + 2 * x.items[0] * x.items[1];

	//return 100 * pow(x.items[1] - x.items[0] * x.items[0], 2) + pow(1 - x.items[0], 2);
}

double funcForSearch(double l, Vector<double> xc, matrix etta, Vector<double> gr) {
	fcount++;
	return func(xc - etta * l * gr);
}


Vector<double> grad(Vector <double> x) {
	Vector<double> t(2);
	t.items[0] = 4 * (x.items[0] - 3) * exp(-pow(x.items[0] - 3, 2) - pow(x.items[1] - 2, 2)) + 
		6 * (x.items[0] - 1) * exp(-pow(x.items[0] - 1, 2) - pow(x.items[1] - 1, 2) / 9);
	t.items[1] = 4 * (x.items[1] - 2) * exp(-pow(x.items[0] - 3, 2) - pow(x.items[1] - 2, 2)) +
		2.0/3.0 * (x.items[1] - 1) * exp(-pow(x.items[0] - 1, 2) - pow(x.items[1] - 1, 2) / 9);

	//t.items[0] = 8 * x.items[0] + 2 * x.items[1];
	//t.items[1] = 2 * x.items[1] + 2 * x.items[0];

	//t.items[0] = 400 * x.items[0] * x.items[0] * x.items[0] - 400 * x.items[0] * x.items[1] + 2 * x.items[0] - 2;
	//t.items[1] = 200 * (x.items[1] - x.items[0] * x.items[0]);
	return t;
}

matrix mull(const Vector<double> & a, const Vector<double> & b) {
	matrix temp(a.size());
	for (int i = 0; i < temp.size; i++) {
		for (int j = 0; j < temp.size; j++) {
			temp.data[i][j] = a.items[i] * b.items[j];
		}
	}
	return temp;
}

void searchMinimumInterval(double & begin, double & end, double x0, double delta, Vector<double> xc, matrix etta, Vector<double> gr)
{
	double f0 = funcForSearch(x0, xc, etta, gr);
	double x1;
	double f1;
	if (f0 < funcForSearch(x0 + delta, xc, etta, gr)) {
		delta = -delta;
	}
	x1 = x0 + delta;
	f1 = funcForSearch(x1, xc, etta, gr);
	while (true) {
		delta *= 2;
		x0 = x1 + delta;
		f0 = funcForSearch(x0, xc, etta, gr);
		if (f1 > f0) {
			x1 = x0;
			f1 = f0;
		}
		else
		{
			x1 -= delta / 2.0;
			break;
		}
	}
	if (x0 > x1) std::swap(x0, x1);
	begin = x0;
	end = x1;
}

double goldenSection(double begin, double end, double eps, Vector<double> xc, matrix etta, Vector<double> gr)
{
	double gold = (3 - sqrt(5.0)) / 2;
	double x1, x2;
	double f1, f2;
	x1 = begin + gold * (end - begin);
	x2 = end - gold * (end - begin);
	f1 = funcForSearch(x1, xc, etta, gr);
	f2 = funcForSearch(x2, xc, etta, gr);
	while (end - begin > eps) {
		if (f1 > f2) {
			begin = x1;
			x1 = x2;
			f1 = f2;
			x2 = end - gold * (end - begin);
			f2 = funcForSearch(x2, xc, etta, gr);
		}
		else
		{
			end = x2;
			x2 = x1;
			f2 = f1;
			x1 = begin + gold * (end - begin);
			f1 = funcForSearch(x1, xc, etta, gr);
		}
	}
	return (begin + end) * 0.5;
}

void broidenMethod(Vector<double>& x, double eps, int maxiter)
{
	std::ofstream out;
	out.open("out11.txt");
	int n = x.size();
	Vector<double> grad1 = grad(x);
	matrix etta(n);
	etta.data[0][0] = 1;
	etta.data[1][1] = 1;
	int iter = 0;
	while (grad1.length() > eps && iter < maxiter) {
		double a, b;
		searchMinimumInterval(a, b, 0, 1e-10, x, etta, grad1);
		double lambda = goldenSection(a, b, 1e-10, x, etta, grad1);
		Vector<double> xk = x - (etta * grad1) * lambda;
		Vector<double> grad2 = grad(xk);
		Vector<double> dgrad = grad2 - grad1;
		Vector<double> dx = xk - x;
		double den = ((dx - etta * dgrad) * dgrad);
		if (abs(den) < EPS) break;
		etta = etta + mull((dx - etta * dgrad), (dx - etta * dgrad)) * (1.0 / den);
		x = xk;
		grad1 = grad2;
		iter++;
		out << iter << '\t' << fcount << '\t' << std::setprecision(14) << x.items[0] << '\t' << x.items[1] << '\t' << func(x) << std::endl;
		if (iter % 1000 == 0) {
			etta.data[0][0] = 1;
			etta.data[0][1] = 0;
			etta.data[1][1] = 1;
			etta.data[1][0] = 0;
			grad1 = grad(x);
		}
	}
}
