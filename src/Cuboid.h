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
    Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
           double m, double distance, double avgVelocityBrownian);

public:
    /**
     * @brief Create particles that belong to the cuboid and insert them into the given vector.
     * @param particles: The vector that particles of the cuboid should be added to.
     */
    void createParticles(std::vector<Particle>& particles) const;
};


