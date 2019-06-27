/**
 * A library for the definitions of the various velocity fields and its gradients.
 */

#include "vecMath3D.hpp"
#include "velocityFields.hpp"

/**
 * The functional definition of the shear field
 * It evalutes the shear field at the given coordinates with the given parameter.
 *
 * @param x, y, z The coordinates of the point
 * @param v0 The scaling factor of the shear field
 * @return The velocity vector at the given point
 */
std::array<double, 3> shearField(double x, double y, double z, double v0) {
	std::array<double, 3> tempReturn = {-sin(M_PI*x)*cos(M_PI*y), cos(M_PI*x)*sin(M_PI*y), 0};
	return v0*tempReturn;
}

/**
 * The functional defintion of the jacobian matrix of the shear field
 * It evaluates the jacobian matrix of the shear field a tthe given coordinates with the given parameter
 *
 * @param x, y, z The coordinates of the point
 * @param v0 The scaling factor of the shear field
 * @return The jacobian matrix at the given point
 */
std::array<std::array<double, 3>, 3> gradShearField(double x, double y, double z, double v0) {
    std::array<std::array<double, 3>, 3> tempReturn;
    tempReturn[0] = {-M_PI*cos(M_PI*x)*cos(M_PI*y), M_PI*sin(M_PI*x)*sin(M_PI*y), 0};
    tempReturn[1] = {-M_PI*sin(M_PI*x)*sin(M_PI*y), M_PI*cos(M_PI*x)*cos(M_PI*y), 0};
    tempReturn[2] = {0, 0, 0};

    tempReturn[0] = v0*tempReturn[0];
    tempReturn[1] = v0*tempReturn[1];
    return tempReturn;
}

/**
 * The functional definition of the navier field
 * It evalutes the navier field at the given coordinates with the given parameters.
 *
 * @param x, y, z The coordinates of the point
 * @param v0 The value of the x-component at the origin
 * @param c1 A parameter of the field
 * @param c2 A parameter of the field
 * @return The velocity vector at the given point
 */
std::array<double, 3> navierField(double x, double y, double z, double v0, double c1, double c2) {
	return { v0 + c1*x + c2*y, -c1*y, 0};
}

/**
 * The functional definition of the navier field
 * It evalutes the navier field at the given coordinates with the given parameters.
 *
 * @param x, y, z The coordinates of the point
 * @param v0 The value of the x-component at the origin
 * @param c1, c2 Parameters of the navier field
 * @return The jacobian matrix at the given point
 */
std::array<std::array<double, 3>, 3> gradNavierField(double x, double y, double z, double v0, double c1, double c2) {
	std::array<std::array<double, 3>, 3> tempReturn;
	tempReturn[0] = {c1, c2, 0};
	tempReturn[1] = {0, -c1, 0};
	tempReturn[2] = {0, 0, 0};
	return tempReturn;
}
