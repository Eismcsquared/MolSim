#include <array>
#include <vector>
#include "body/Particle.h"

/**
 * @brief This class represents a cell in the simulation domain that is used in the linked cell algorithm.
 */
class Cell {
private:

   /**
    * The position of the left (negative x-direction) lower (negative y-direction) back-side (negative z-direction) corner
    */
    std::array<double, 3> position;

   /**
    * The size of the cell, given as width (length in x-direction), height (length in y-direction) and depth (length in z-direction).
    */
    std::array<double, 3> size;

   /**
    * The indices of particles contained in a cell.
    */
    std::vector<int> particleIndices;

public:
   /**
    * Constructor.
    * @param position The position of the cell.
    * @param size The size of the cell.
    */
    Cell(std::array<double, 3> position, std::array<double, 3> size);
   /**
    * The getter for the position of a cell.
    * @return The position of the cell.
    */
    const std::array<double, 3> getPosition() const;
   /**
    * The getter for the size of a cell.
    * @return The size of the cell.
    */
    const std::array<double, 3> getSize() const;
   /**
    * The getter for the particles' indices in a cell.
    * @return The particles' indices in the cell as a vector.
    */
    const std::vector<int> &getParticleIndices() const;
   /**
    * Test whether a given position belongs to the cell.
    * @param position The position that should be tested.
    * @return True if the given position is in the cell.
    */
    bool contains(std::array<double, 3> pos);
   /**
    * Add a index to the cell.
    * @param index The index to be added into the cell.
    */
    void addIndex(int index);
   /**
    * Remove a index to the cell.
    * @param index The index to be removed from the cell.
    */
    void removeIndex(int index);

};

