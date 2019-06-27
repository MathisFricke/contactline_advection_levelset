#ifndef FIELDS_HPP_
#define FIELDS_HPP_
#include "vecMath3D.hpp"

std::array<double, 3> shearField(double x, double y, double z, double v0);
std::array<std::array<double, 3>, 3> gradShearField(double x, double y, double z, double v0);
std::array<double, 3> navierField(double x, double y, double z, double v0, double c1, double c2);
std::array<std::array<double, 3>, 3> gradNavierField(double x, double y, double z, double v0, double c1, double c2);

#endif /* FIELDS_HPP_ */
