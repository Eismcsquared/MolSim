#include <array>
#include <vector>
#include "Particle.h"
#include "Cluster.h"
/**
 * @brief This class represents a cuboid consisting of particles.
 */
class Cuboid : public Cluster {
private:
    /**
     * The number of particles in x, y and z direction.
     */
    std::array<unsigned int, 3> N;

public:
    Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
           double m, double distance, double avgVelocityBrownian, int dimension);

    /**
     * Create particles that belong to the cuboid and insert them into the given vector.
     * @param particles: The vector that particles of the cuboid should be added to.
     */
    void createParticles(std::vector<Particle>& particles) const;
};


