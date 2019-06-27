#ifndef CLASS_LEVELSET
#define CLASS_LEVELSET

#include <sstream>    // Stringstream
#include <fstream>    // Filestream

#include "Field.hpp"
#include "VelocityField.hpp"
#include "vecMath3D.hpp"

class LevelSet : Field<double> {
private:
    //The width of a cell in each direction
    double dx, dy, dz;

    // A pointer to the VelocityField acting on the LevelSet field.
    VelocityField *field;

    // Decides which contact point to track. Only applicable in 2D.
    std::string trackedCP;

public:
    LevelSet(int numX, int numY, int numZ, double dx, double dy, double dz, VelocityField *field, std::string trackedCP);

    array<double, 3> getInitCP(array<double, 3> expcp, double epsilon);
    array<double, 3> getContactPointExplicitEuler(double dt, int timestep, array<double, 3> initCP);
    array<double, 3> getContactPointLinearField(double t, double c1, double x0, double v0);
    array<int, 3> getContactPointIndices(array<double, 3> point);
    array<double, 3> getContactPoint(array<int, 3> indices);
    double getReferenceAngleExplicitEuler(double dt, int timestep, array<double, 3> n_sigma_init, array<double, 3> CP_init);
    double getReferenceAngleLinearField(double t, double c1, double c2, double theta0);
    double getContactAngle(array<int, 3> cell);
    double sumLevelSet();
    void writeToFile(double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *xmfFile);
    void initDroplet(array<double, 3> center, double radius);
    void calculateNextTimestep(double dt, int timestep);
};

#endif
