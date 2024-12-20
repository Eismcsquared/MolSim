#pragma once

#include <vector>
#include "Particle.h"
#include "utils/MaxwellBoltzmannDistribution.h"
/**
 * @brief The class represents a cluster of particles
 */
class Cluster {
protected:
    /**
     * The position of the cluster:
     * For cuboid: the left (negative x-direction) lower (negative y-direction) back-side (negative z-direction) corner of the cuboid.
     * For sphere: the center of the sphere
     */
    std::array<double, 3> x;
    /**
     * The velocity of the cluster.
     */
    std::array<double, 3> v;

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
    /**
     * The dimension of the cluster, either 2 or 3.
     */
     int dimension;

public:

    /**
     * Constructor that initializes all the attributes.
     */
    Cluster(const std::array<double, 3> &x, const std::array<double, 3> &v,
           double m, double distance, double avgVelocityBrownian, int dimension);

    /**
     * Create particles that belong to the cuboid and insert them into the given vector.
     * @param particles: The vector that particles of the cuboid should be added to.
     */
    virtual void createParticles(std::vector<Particle>& particles) const = 0;
};