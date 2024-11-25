#include <array>
#include <vector>
#include "body/Particle.h"
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
    * The particles contained in a cell.
    */
    std::vector<Particle> particles;

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
    const std::array<double, 3> getPosition();
   /**
    * The getter for the size of a cell.
    * @return The size of the cell.
    */
    const std::array<double, 3> getSize();
   /**
    * The getter for the particles in a cell.
    * @return The particles in the cell as a vector.
    */
    const std::vector<Particle> getParticles();
   /**
    * Test whether a given position belongs to the cell.
    * @param position The position that should be tested.
    * @return True if the given position is in the cell.
    */
    bool contains(std::array<double, 3> pos);
   /**
    * Add a particle to the cell. Nothing happens if the position of the particle does not belong to the cell.
    * @param particle The particle to be added into the cell.
    */
    void addParticle(Particle& particle);
   /**
    * Remove a particle to the cell.
    * @param particle The particle to be removed from the cell.
    */
    void removeParticle(Particle& particle);

};

