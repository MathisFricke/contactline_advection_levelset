/**
 * @class Field
 * The template for the abstract field class.
 *
 * Save the attribute T for each point on the grid
 */
#ifndef CLASS_FIELD
#define CLASS_FIELD
#include <vector>

template <class T>
class Field {
    private:
    //Used to store the data of the field
    std::vector<T> data;

protected:
    /// The number of cells in each direction
    ///@{
    int numX, numY, numZ;
    ///@}
public:
    Field(int numX, int numY, int numZ);
    T& at(int x, int y, int z);
    const T& at(int x, int y, int z) const;

};
#include "Field.cpp"

#endif
