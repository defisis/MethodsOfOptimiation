#pragma once
#include <iostream>
#include <iomanip>

template <typename T>
struct Vector {
private:
	int size_;

public:
	T * items = nullptr;

	Vector();
	Vector(int size);
	Vector(int size, T _item);
	Vector(const Vector <T> & v);
	~Vector();

	double length() const;
	int size() const;
	void resize(int size);
	const Vector<T> operator+(const Vector<T> &op) const;
	const Vector<T> operator-(const Vector<T> &op) const;
	const Vector<T> operator*(const T &op) const;
	const double operator*(const Vector<T> &op) const;
	Vector<T>& operator=(const Vector<T>& op);

	friend const Vector<T> operator*(const double &op1, const Vector<T> &op2)
	{
		return op2 * op1;
	}
	friend Vector<T>& operator*=(Vector<T> &op1, const T &op2)
	{
		op1 = op1 * op2;
		return op1;
	}
	friend Vector<T>& operator+=(Vector<T> &op1, const Vector<T> &op2)
	{
		op1 = op1 + op2;
		return op1;
	}
	friend Vector<T>& operator-=(Vector<T> &op1, const Vector<T> &op2)
	{
		op1 = op1 - op2;
		return op1;
	}
	friend std::istream& operator >> (std::istream &in, const Vector<T> &obj)
	{
		for (int i = 0; i < obj.size(); i++) {
			in >> obj.items[i];
		}
		return in;
	}
	friend std::ostream& operator << (std::ostream &out, const Vector<T> &obj)
	{
		for (int i = 0; i < obj.size(); i++) {
			out << std::setprecision(16) << obj.items[i] << '\n';
		}
		return out;
	}
};

template <typename T>
Vector<T>::Vector()
{
}

template<typename T>
inline Vector<T>::Vector(int size) : size_{ size }
{
	items = new T[size_];
}

template<typename T>
inline Vector<T>::Vector(int size, T _item) : size_{ size }
{
	items = new T[size_];
	for (int i = 0; i < size_; i++) {
		items[i] = _item;
	}
}

template<typename T>
inline Vector<T>::Vector(const Vector<T> & v) : size_{ v.size_ }
{
	if (!items)
		items = new T[size_];
	for (int i = 0; i < size_; i++) {
		items[i] = v.items[i];
	}
}

template<typename T>
inline Vector<T>::~Vector()
{
	if (items) {
		delete [] items;
		items = nullptr;
	}
}

template<typename T>
inline double Vector<T>::length() const
{
	double sum = 0;
	for (int i = 0; i < size_; i++) {
		sum += items[i] * items[i];
	}
	return sqrt(sum);
}

template<typename T>
inline int Vector<T>::size() const
{
	return size_;
}

template<typename T>
inline void Vector<T>::resize(int size)
{
	if (items) {
		delete[] items;
	}
	items = new T[size];
	size_ = size;
}


template <typename T>
const Vector<T> Vector<T>::operator+(const Vector<T> &op) const
{
	Vector _temp(size_);
	for (int i = 0; i < size_; i++) {
		_temp.items[i] = items[i] + op.items[i];
	}
	return _temp;
}

template <typename T>
const Vector<T> Vector<T>::operator-(const Vector<T> &op) const
{
	Vector _temp(size_);
	for (int i = 0; i < size_; i++) {
		_temp.items[i] = items[i] - op.items[i];
	}
	return _temp;
}

template <typename T>
const Vector<T> Vector<T>::operator*(const T & op) const
{
	Vector _temp(size_);
	for (int i = 0; i < size_; i++) {
		_temp.items[i] = items[i] * op;
	}
	return _temp;
}

template<typename T>
const double Vector<T>::operator*(const Vector<T>& op) const
{
	double _temp_sum = 0;
	for (int i = 0; i < size_; i++) {
		_temp_sum += items[i] * op.items[i];
	}
	return _temp_sum;
}

template<typename T>
inline Vector<T>& Vector<T>::operator=(const Vector<T>& op)
{
	if (this == &op)
		return *this;
	else {
		for (int i = 0; i < size_; i++) {
			items[i] = op.items[i];
		}
	}
	return *this;
}

