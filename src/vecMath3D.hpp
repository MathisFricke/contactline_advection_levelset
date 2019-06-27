#ifndef VEC_MATH_3D
#define VEC_MATH_3D

#include <array>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

using std::array;

array<double, 3> operator+ (array<double, 3> vecA, array<double, 3> vecB);

array<int, 3> operator+ (array<int, 3> vecA, array<int, 3> vecB);

array<double, 3> operator- (array<double, 3> vecA, array<double, 3> vecB);

array<double, 3> operator- (array<double, 3> vecA, array<int, 3> vecB);

array<double, 3> operator*(double a, const array<double, 3>& vec);

array<double, 3> operator*(const array<double, 3>& vec, double a);

array<double, 3> operator*(const array<int, 3>& vec, double a);

array<array<double, 3>, 3> operator* (double a, const array<array<double,3>, 3>& matrix);

double operator* (const array<double, 3>& vecA, const array<double, 3>& vecB);

array<double, 3> operator* (const array<array<double,3>, 3>& matrix, const array<double, 3> vec);

array<double, 3> operator/ (const array<double, 3>& vec, double a);

double abs (const array<double, 3> vec);

array<array<double, 3>, 3> transpose(const array<array<double, 3>, 3> matrix);

#endif
