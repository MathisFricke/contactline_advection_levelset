#ifndef CLASS_VELOCITYFIELD
#define CLASS_VELOCITYFIELD
#include <array>

using std::array;

class VelocityField {
private:
    //@{
    /** The space the velocity field is defined on */
	double xmin, xmax, ymin, ymax, zmin, zmax;
	//@}

	//@{
	/** The width of a cell in each direction */
	double dx, dy, dz;
	//@}

	//@{
	/** The parameters of the velocity field */
	double v0, c1, c2, tau;
	//@}

	/// The kind of velocity field, either navier, navier with cosine modulation or shear.
	std::string name;
	///The maximum absolute value of the field on the space it is defined on.
	double maxAbsoluteValue;
public:
	VelocityField(std::string name, double v0, double c1, double c2, double tau,
			double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double dx, double dy, double dz);
	array<double, 3> at(double t, double x, double y, double z);
	array<array<double, 3>, 3> gradAt(double t, double x, double y, double z);
	void writeToFile(double t);
	double getC1();
	double getTau();
	std::string getName();
	double getMaxNormValue();
};

#endif
