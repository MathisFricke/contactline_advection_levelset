/**
 * @class LevelSet
 * The class of the Level set field.
 *
 * This class contains most of the methods used in the program, such as calculating the reference contact angle,
 * calculating the contact angle from the simulation data, tracking the contact point and evolving the field with time.
 */

#include "LevelSet.hpp"

/**
 * The constructor.
 * @param numX The number of cells in x direction
 * @param dx Width of one cell in x direction
 * @param field Pointer to the velocity field object acting on the levelset field
 * @param trackedCP Which contact point to track. Only applicable in 2D.
 */
LevelSet::LevelSet(int numX, int numY, int numZ, double dx, double dy, double dz, VelocityField *field, std::string trackedCP) : Field<double>(numX, numY, numZ) {
		this->dx = dx;
		this->dy = dy;
		this->dz = dz;
		this->field = field;
		this->trackedCP = trackedCP;
}

/**
 * Calculate the initial contact point from a given expected contact point.
 * Using the expected contact point, iterate over all points of the field, looking for the closest point to expcp which
 * has an absolute Level set value of less than epsilon and lies of on the interface.
 *
 * @param expcp The expected contact point
 * @param epsilon Maximum absolute Level set value to consider a point "on the interface"
 * @return The contact point closest to the expected contact point
 *
 */
array<double, 3> LevelSet::getInitCP(array<double, 3> expcp, double epsilon) {
array<double, 3> candidate = {0, 0, 0};
    for (int x = 0; x < this->numX; x++)
        for (int y = 0; y < this->numY; y++)
            for (int z = 0; z < this->numZ; z++) {
                array<double, 3> other = {x*dx, y*dy, z*dz};
                if (abs(candidate - expcp) > abs(other - expcp) && std::abs(this->at(x, y, z)) < epsilon && y == 0) {
                    candidate = other;
		}
	}

    return candidate;
}

/**
 * Calculate the new position of the contact point at a given time.
 * Using the explicit euler method, calculate the position of the contact point initCP at time timestep*dt.
 *
 * @param dt The length of a single timestep
 * @param timestep The index of the timestep
 * @param initCP The position of the contact point at timestep 0
 * @return The theoretical position of the contact point after a given time
 */
array<double, 3> LevelSet::getContactPointExplicitEuler(double dt, int timestep, array<double, 3> initCP) {
    array<double, 3> &temp = initCP;
    for (int i = 0; i < timestep; i++)
    	temp = temp + dt*field->at(timestep*dt, temp[0], temp[1], temp[2]);

    return temp;
}

/**
 * Calculate the contact point position for the navier field using the analytic solution.
 *
 * @param t The time at which to calculate the contact point position
 * @param c1 A parameter of the navier field
 * @param x0 The initial position of the contact point
 * @param v0 A parameter of the navier field
 * @return The position of the contact point in arbitrary units
 */
array<double, 3> LevelSet::getContactPointLinearField(double t, double c1, double x0, double v0) {
	if (field->getName() == "navierField") {
		return {x0 * exp(c1 * t) + v0/c1 * (exp(c1 * t) - 1), 0, 0};
	} else if (field->getName() == "timeDependentNavierField") {
		double tau = field->getTau();
		return {x0 *exp(c1/M_PI * tau*sin(M_PI*t/tau)) + v0/c1 * (exp(c1/M_PI * tau*sin(M_PI*t/tau)) -1), 0, 0};
	} else {
		throw std::invalid_argument("Please choose either navierField or timeDependentNavierField if analyzing linear fields.");
	}
}

/**
 * Return the indices of the contact point in 2D or the coordinates most closely matching the given parameter.
 * This function works differently for 2D and 3D. If the simulation is 2D, it ignores the input parameter and returns
 * indices of the contact point by checking where the sign of the LevelSet field changes.
 * In 3D, it returns the coordinates of the point on the grid which is closest to the parameter.
 *
 * @param point A point in the space of the Level set field.
 * @return 2D: Indices of the contact point. 3D: Indices of the point closest to point.
 */
array<int, 3> LevelSet::getContactPointIndices(array<double, 3> point) {
    if (this->numZ == 1){
        if (trackedCP == "left") {
            double initSign = this->at(0, 0, 0)/std::abs(this->at(0, 0, 0));
            for (int i = 0; i < numX; i++)
                if (this->at(i, 0, 0)*initSign < 0) {
                	return {i, 0, 0};
                }
        } else {
            double initSign = this->at(numX - 1, 0, 0)/std::abs(this->at(numX - 1, 0, 0));
            for (int i = numX - 1; i >= 0; i--)
                if (this->at(i, 0, 0)*initSign < 0) {
                	return {i, 0, 0};
				}
		}
    } else {
        //Find cell corresponding to this point
        array<int, 3> cell = {0, 0, 0};
        for (int x = 0; x < this->numX; x++)
            for (int y = 0; y < this->numY; y++)
                for (int z = 0; z < this->numZ; z++) {
                    array<int, 3> other = {x, y, z};
                    if (abs(point - cell*dx) > abs(point - other*dx))
                        cell = other;
                }
        return cell;
    }
    return array<int, 3> {0,0,0};
}

/**
 * Get contact point coordinates from contact point indices.
 *
 * Uses a convex combination of neighboring LevelSet values depending on whether the right or the
 * left contact point is being tracked.
 *
 * @param indices The indices of the cell where the first change of sign occured during iteration
 * inside LevelSet::getContactPointIndices
 */
array<double, 3> LevelSet::getContactPoint(array<int, 3> indices) {
	int i = indices[0];
	if (trackedCP == "left") {
		double alpha = this->at(i,0,0)-this->at(i -1,0,0);

		if(std::abs(alpha) < 1E-12){
			throw std::runtime_error("Difference of LevelSet values at contact point too small for convex combination");
		}

		alpha = this->at(i,0,0)/(alpha);
		return {(1 - alpha)*i*dx + alpha*(i - 1)*dx, 0, 0};

	} else if (trackedCP == "right") {
		double alpha = this->at(i + 1, 0, 0) - this->at(i, 0, 0);
		if(std::abs(alpha) < 1E-12) {
		   throw std::runtime_error("Difference of LevelSet values at contact point too small for convex combination");
		}

		alpha = this->at(i + 1, 0, 0) / alpha;
		return {alpha*i*dx + (1 - alpha)*(i+1)*dx, 0, 0};

	} else {
		throw std::runtime_error("Variable \"trackedContactPoint\" is not set to left or right");
	}
}

/**
 * Calculate the contact angle.
 * This function uses finite differences to calculate the normal vector of the Level set field at cell and
 * thus calculate the contact angle.
 *
 * @param cell The indices of the contact point.
 * @return The contact angle in degrees
 */
double LevelSet::getContactAngle(array<int, 3> cell) {

	double normalX = 0, normalY = 0, normalZ = 0;

	if (trackedCP == "left") {
		 // find root of phi, alpha: coefficient for convex combination
		 double alpha = this->at(cell[0],0,0)-this->at(cell[0]-1,0,0);

		 if(std::abs(alpha) < 1E-12){
		   throw std::runtime_error("Difference of LevelSet values at contact point too small for convex combination");
		 }

		 alpha = this->at(cell[0],0,0)/(alpha);

		//Calculate angle at this cell with finite differences
		normalX = alpha*(this->at(cell[0], cell[1], cell[2]) - this->at(cell[0]-2, cell[1], cell[2]))/(2*dx)
		+ (1-alpha)*(this->at(cell[0]+1, cell[1], cell[2]) - this->at(cell[0]-1, cell[1], cell[2]))/(2*dx);
		normalY = alpha*(-this->at(cell[0]-1, cell[1] + 2, cell[2])
						  + 4.0*this->at(cell[0]-1, cell[1] + 1, cell[2])
						  - 3.0*this->at(cell[0]-1, cell[1], cell[2]))/(2*dy)
						  + (1-alpha)*(-this->at(cell[0], cell[1] + 2, cell[2])
						  + 4.0*this->at(cell[0], cell[1] + 1, cell[2])
						  - 3.0*this->at(cell[0], cell[1], cell[2]))/(2*dy); // second order difference quotient

	} else if (trackedCP == "right") {
		double alpha = this->at(cell[0] + 1, 0, 0) - this->at(cell[0], 0, 0);
        if(std::abs(alpha) < 1E-12) {
		   throw std::runtime_error("Difference of LevelSet values at contact point too small for convex combination");
        }

        alpha = this->at(cell[0] + 1, 0, 0) / alpha;

        normalX = alpha*(this->at(cell[0] + 1, 0, 0) - this->at(cell[0] - 1, 0, 0))/(2*dx)
        		+ (1 - alpha)*(this->at(cell[0] + 2, 0, 0) - this->at(cell[0], 0, 0))/(2*dx);
        normalY = alpha*(-this->at(cell[0], cell[1] + 2, cell[2])
						  + 4.0*this->at(cell[0], cell[1] + 1, cell[2])
						  - 3.0*this->at(cell[0], cell[1], cell[2]))/(2*dy)
						  + (1 - alpha)*(-this->at(cell[0] + 1, cell[1] + 2, cell[2])
						  + 4.0*this->at(cell[0] + 1, cell[1] + 1, cell[2])
						  - 3.0*this->at(cell[0] + 1, cell[1], cell[2]))/(2*dy);
	}
    if (this->numZ > 1)
    	normalZ = (this->at(cell[0], cell[1], cell[2]+1) - this->at(cell[0]-1, cell[1], cell[2]-1))/(2*dz);
    else
    	normalZ = 0;
    array<double ,3> normal = {normalX, normalY, normalZ};
    normal = normal/abs(normal);
    return acos(normal[1]);
}

/**
 * Calculate the contact angle at a given time using the explicit Euler method.
 *
 * @param dt The length of a timestep
 * @param timestep The index of the timestep
 * @param n_sigma_init The initial normal vector of the interface
 * @param CP The current position of the contact point
 * @return The reference contact angle in degrees
 */
double LevelSet::getReferenceAngleExplicitEuler(double dt, int timestep, array<double, 3> n_sigma_init, array<double, 3> CP_init) {
	array<double, 3> &n_sigma = n_sigma_init;
	array<double, 3> deriv = {0, 0, 0};
	double t = dt*timestep;
	for (int i = 0; i < timestep; i++) {
	    array<double, 3> CP = getContactPointExplicitEuler(dt, i, CP_init);
		deriv = -1*transpose(field->gradAt(t, CP[0], CP[1], CP[2]))*n_sigma + ((field->gradAt(t, CP[0], CP[1], CP[2])*n_sigma)*n_sigma)*n_sigma;
		n_sigma = deriv*dt + n_sigma;
		n_sigma = n_sigma/abs(n_sigma);
	}
	// Normalize n_sigma. In theory this should not be necessary since the vector should stay normalized, although it may reduce numerical inaccuracies
	return acos(n_sigma[1]);
}

/**
 * Calculate the contact angle at a given time for the navier field using the analytic solution.
 *
 * @param t The time at which to calculate the contact angle
 * @param c1 A parameter of the navier field
 * @param c2 A parameter of the navier field
 * @param theta0 The contact angle at time t == 0
 * @return The reference contact angle in degrees
 */
double LevelSet::getReferenceAngleLinearField(double t, double c1, double c2, double theta0) {
	if (field->getName() == "navierField")
	    if (trackedCP == "left") {
	        return M_PI/2 + atan(-1/tan(theta0) * exp(2*c1*t) - c2 * (exp(2*c1*t) - 1)/(2*c1));
	    } else {
	        return M_PI/2 + atan(-1/tan(theta0) * exp(2*c1*t) + c2 * (exp(2*c1*t) - 1)/(2*c1));
	    }
	else if (field->getName() == "timeDependentNavierField") {
		double tau = field->getTau();
		if (trackedCP == "left") {
		    return M_PI/2 + atan(-1/tan(theta0) * exp(2*c1*tau/M_PI*sin(M_PI*t/tau)) - c2 * (exp(2*c1*tau/M_PI*sin(M_PI*t/tau)) - 1)/(2*c1));
		} else {
		    return M_PI/2 + atan(-1/tan(theta0) * exp(2*c1*tau/M_PI*sin(M_PI*t/tau)) + c2 * (exp(2*c1*tau/M_PI*sin(M_PI*t/tau)) - 1)/(2*c1));
		}
	}
	else
		throw std::invalid_argument("Please choose either navierField or timeDependentNavierField if analyzing linear fields.");
}

/**
 * Write the files necessary to visualize the field and its velocity field in Paraview.
 *
 * Write the XMF file to disk, as well as the grid of the Level set field. This is done only once.
 * At each point in time, write the Level set field to disk and call VelocityField::writeToFile to write its values as well.
 *
 * @param dt The width of a timestep
 * @param timestep The index of the timestep
 * @param total_timesteps The total number of timesteps
 * @param total_writesteps The total number of steps that will be written to disk
 * @param xmfFile A pointer to the XMF file.
 */
void LevelSet::writeToFile(double dt, int timestep, int total_timesteps, int total_writesteps, std::ofstream *xmfFile) {
    int Npoints = numX*numY*numZ;
    double *pointCoordinates = new double[Npoints*3];
    double *pointPhiValues = new double [Npoints];

    int index = 0;
    for (int k = 0; k < numZ; k++)
        for (int j = 0; j < numY; j++)
            for (int i = 0; i < numX; i++) {
                if (timestep == 0) {
                    pointCoordinates[index] = i*dx;
                    pointCoordinates[index +1] = j*dy;
                    pointCoordinates[index +2] = k*dz;
                    index += 3;
                }
                pointPhiValues[i + j*numX + k*numX*numY] = this->at(i, j, k);
	    }

    // If it is the first iteration, create coordinate file
    if (timestep == 0) {
        *xmfFile << "<?xml version=\"1.0\" ?>\n"
                 << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
                 << "<!ENTITY Npoints \"" + std::to_string(Npoints) + "\">\n"
                 << "]>\n"
                 << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n"
                 << "<Domain>\n"
                 << "<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
                 << "<Time TimeType=\"List\">\n"
                 << "<DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\""+ std::to_string(total_writesteps)+"\">\n";
        for (int i = 0; i < total_writesteps; i++) {
            *xmfFile << std::to_string(i*((double)total_timesteps/total_writesteps)*dt) << " ";
        }
        *xmfFile <<"</DataItem>\n" << "</Time>\n";

        //Write field coordinates into binary file
        FILE *fieldFile;
        fieldFile = fopen("data/field.bin", "wb");
        fwrite(pointCoordinates, sizeof(double), Npoints*3, fieldFile);
        fclose(fieldFile);
    }
    FILE *PhiFile;
    std::string filename = "data/Phi_t="+ std::to_string(timestep*dt)+".bin";
    PhiFile = fopen(filename.data(), "wb");
	fwrite(pointPhiValues, sizeof(double), Npoints, PhiFile);
	fclose(PhiFile);



	*xmfFile << "<Grid>\n"
			 << "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\""+std::to_string(Npoints) +"\"/>\n"
			 << "<Geometry GeometryType=\"XYZ\"> \n"
			 << "<DataItem Name=\"points\" Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Endian=\"Little\" Dimensions=\"&Npoints; 3\">\n"
			 << "field.bin\n"
			 << "</DataItem>\n</Geometry>\n"
			 << "<Attribute Name =\"lvlset\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
			 << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\"  Endian=\"Little\" Dimensions=\"&Npoints;\">\n"
			 << "Phi_t=" + std::to_string(timestep*dt) +".bin\n"
			 << "</DataItem></Attribute>\n"
			 << "<Attribute Name =\"VelocityField\" AttributeType=\"Vector\" Center=\"Cell\">\n"
			 << "<DataItem Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Endian=\"Little\" Dimensions=\"&Npoints; 3\">\n"
			 << "Vel_t=" + std::to_string(timestep*dt) +".bin\n"
			 << "</DataItem></Attribute>\n"
			 <<"</Grid>\n";


    delete[] pointCoordinates;
    delete[] pointPhiValues;

    //Write velocity field
	field->writeToFile(timestep*dt);
}

/**
 * Calculate the sum of the Level set field.
 * @return The value of the sum
 */
double LevelSet::sumLevelSet() {
    double temp = 0;
    for (int x = 0; x < this->numX; x++)
        for (int y = 0; y < this->numY; y++)
            for (int z = 0; z < this->numZ; z++)
                temp = temp + this->at(x, y, z);

    return temp;
}

/**
 * Initialize the droplet.
 * The droplet is initialized on the Level set field by setting the space occupied by the fluid to be negative,
 * the space occupied by the gas to be positive and the zero set to be the interface. More specifically, each
 * point on the grid is set to the value
 * \f[ \Phi = (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 \f]
 * Where \f$ x_0, y_0, z_0 \f$ are the coordinates of the droplet center and R is its radius.
 *
 * @param center The center of the droplet
 * @param radius The radius of the droplet
 */
void LevelSet::initDroplet(array<double, 3> center, double radius) {
    for (int x = 0; x < this->numX; x++)
        for (int y = 0; y < this->numY; y++)
            for (int z = 0; z < this->numZ; z++)
                this->at(x, y, z) = pow(x*dx - center[0], 2) + pow(y*dy - center[1], 2) + pow(z*dz - center[2], 2) - pow(radius, 2);
}

/**
 * Calculate the next timestep at a given time and evolves the field.
 *
 * For the flux calculation, the upwind method is used to ensure stability. Furthermore we impose
 * the Neumann boundary condition onto the field
 *
 * @param dt The width of the timestep by which to evolve the field
 * @param timestep The index of the timestep
 */
void LevelSet::calculateNextTimestep(double dt, int timestep) {
    LevelSet tempPhi(*this);

    const array<double, 3> upNormal = {0, 1, 0};
    const array<double, 3> downNormal = {0,-1, 0};
    const array<double, 3> leftNormal = {-1, 0, 0};
    const array<double, 3> rightNormal = {1, 0, 0};
    const array<double, 3> frontNormal = {0, 0, 1};
    const array<double, 3> backNormal = {0, 0, -1};

    // loop over all cells
    for (int x = 0; x < numX; x++) {
        for (int y = 0; y < numY; y++) {
            for (int z = 0; z < numZ; z++) {

                //Calculate the flux of phi over the cell faces
                double flux = 0;
                double sp = 0;

                for (int dir = 0; dir < 6; dir++) {

                    switch(dir) {
                    case 0:
                        sp = field->at(timestep*dt, (x+1/2)*dx, (y+1)*dy, (z+1/2)*dz) * upNormal;
                        if (y == numY - 1) {
                            flux += sp*tempPhi.at(x, y, z)*dx*dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y+1, z))*dx*dz;
                        }
                        break;

                    case 1:
                        sp = field->at(timestep*dt, (x+1/2)*dx, y*dy, (z+1/2)*dz) * downNormal;
                        if (y == 0 ) {
                            flux += sp*tempPhi.at(x, y, z)*dx*dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x, y-1, z))*dx*dz;
                        }
                        break;

                    case 2:
                        sp = field->at(timestep*dt, x*dx, (y+1/2)*dy, (z+1/2)*dz) * leftNormal;
                        if (x == 0) {
                            flux += sp*tempPhi.at(x, y, z)*dy*dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x-1, y, z))*dy*dz;
                        }
                        break;

                    case 3:
                        sp = field->at(timestep*dt, (x+1)*dx, (y+1/2)*dy, (z+1/2)*dz) * rightNormal;
                        if (x == numX - 1) {
                            flux += sp*tempPhi.at(x, y, z)*dy*dz;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z) + fmin(sp,0.0)*tempPhi.at(x+1, y, z))*dy*dz;
                        }
                        break;

                    case 4:
                        if (numZ == 1)
                            break; // only relevant for 3D
                        sp = field->at(timestep*dt, (x+1/2)*dx, (y+1/2)*dy, (z+1)*dz) * frontNormal;
                        if (z == numZ - 1) {
                            flux += sp*tempPhi.at(x, y, z)*dx*dy;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z+1) + fmin(sp,0.0)*tempPhi.at(x, y, z))*dx*dy;
                        }
                        break;

                    case 5:
                        if (numZ == 1)
                            break; // only relevant for 3D
                        sp = field->at(timestep*dt, (x+1/2)*dx, (y+1/2)*dy, z*dz) * backNormal;
                        if (z == 0) {
                            flux += sp*tempPhi.at(x, y, z)*dx*dy;
                        }
                        else{
                            flux += (fmax(sp,0.0)*tempPhi.at(x, y, z-1) + fmin(sp,0.0)*tempPhi.at(x, y, z))*dx*dy;
                        }
                        break;
                    }
                }
                this->at(x, y, z) = this->at(x, y, z) - dt/(dx*dy*dz)*flux;
            }
        }
    }
}
