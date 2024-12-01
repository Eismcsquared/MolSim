#include <array>
#include "Cluster.h"
/**
 * @brief This class represents a sphere consisting of particles
 */
class Sphere : public Cluster {
    /**
     * The radius of the sphere in terms of number of particles along the radius.
     */
    int radius;

public:
    Sphere(const std::array<double, 3> &x, const std::array<double, 3> &v, int radius,
           double m, double distance, double avgVelocityBrownian, int dimension);

    /**
     * Create particles that belong to the sphere and insert them into the given vector.
     * @param particles: The vector that particles of the cuboid should be added to.
     */
    void createParticles(std::vector<Particle>& particles) const;
};
