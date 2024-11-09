#include <array>
#include <vector>
#include "Particle.h"

class Cuboid {
private:
    /**
     * the lower, left, front corner (i.e. the corner with least coordinates) of the cuboid.
     */
    std::array<double, 3> x;
    /**
     * The velocity of the cuboid.
     */
    std::array<double, 3> v;
    /**
     * The number of particles in x, y and z direction.
     */
    std::array<unsigned int, 3> N;
    /**
     * The mass of each particle.
     */
    double m;
    /**
     * The mesh width of the grid.
     */
    double distance;
    /**
     * The average velocity of the Brownian motion.
     */
    double avgVelocityBrownian;

public:
    std::vector<Particle> createParticles();
};

