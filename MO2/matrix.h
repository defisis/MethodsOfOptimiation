#pragma once
#include "Vector.h"
struct matrix {
	int size = 0;
	double **data = nullptr;
	matrix(int n) {
		size = n;
		data = new double*[n];
		for (int i = 0; i < n; i++) {
			data[i] = new double[n];
			for (int j = 0; j < n; j++) {
				data[i][j] = 0;
			}
		}
	}
	matrix(const matrix & m) {
		size = m.size;
		if (!data) {
			data = new double*[size];
			for (int i = 0; i < size; i++) {
				data[i] = new double[size];
			}
		}
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				data[i][j] = m.data[i][j];
			}
		}
	}
	~matrix() {
		for (int i = 0; i < size; i++) {
			delete[] data[i];
		}
		delete[] data;
		data = nullptr;
	}
	const Vector<double> operator*(const Vector<double> &op) const;
	const matrix operator*(const double &op) const;
	const matrix operator+(const matrix& op) const; 
	matrix& operator=(const matrix& op);
	

	friend const Vector<double> operator*(const Vector<double> &op1, const matrix &op2)
	{
		Vector<double> temp(op1.size());
		for (int i = 0; i < temp.size(); i++) {
			for (int j = 0; j < temp.size(); j++) {
				temp.items[i] = op1.items[j] * op2.data[j][i];
			}
		}
		return temp;
	}
};


inline const Vector<double> matrix::operator*(const Vector<double>& op) const
{
	Vector<double> temp(op.size());
	for (int i = 0; i < size; i++) {
		temp.items[i] = 0;
		for (int j = 0; j < size; j++) {
			temp.items[i] += data[i][j] * op.items[j];
		}
	}
	return temp;
}

inline const matrix matrix::operator*(const double & op) const
{
	matrix temp(size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			temp.data[i][j] = data[i][j] * op;
		}
	}
	return temp;
}

inline const matrix matrix::operator+(const matrix & op) const
{
	matrix temp(size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			temp.data[i][j] = data[i][j] + op.data[i][j];
		}
	}
	return temp;
}

inline matrix & matrix::operator=(const matrix & op)
{
	if (this != &op) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				data[i][j] = op.data[i][j];
			}
		}
	}
	return *this;
}

