#pragma once

#include <array>
#include <vector>
#include "Particle.h"
#include "utils/ArrayUtils.h"
#include "Cluster.h"
/**
 * @brief This class represents a cuboid consisting of particles.
 */
class Cuboid : public Cluster {
protected:
    /**
     * The number of particles in x, y and z direction.
     */
    std::array<unsigned int, 3> N;

    /**
     * Marks the cuboid as stationary (i.e. wall)
     */
    bool stationary;

public:
    /**
     * Constructor. The parameters epsilon and sigma are set to 5 and 1 by default, respectively.
     * @param x The left-down-back side corner (the corner with the least coordinates) of the cuboid.
     * @param v The mean velocity of the cuboid.
     * @param n The number of particles in each direction.
     * @param m The mass of each particle.
     * @param distance The distance between neighbouring particles.
     * @param avgVelocityBrownian The average velocity of the Brownian motion.
     * @param dimension The dimension of the Brownian motion.
     */
    Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
           double m, double distance, double avgVelocityBrownian, int dimension);

    /**
     * Constructor.
     * @param x The left-down-back side corner (the corner with the least coordinates) of the cuboid.
     * @param v The mean velocity of the cuboid.
     * @param n The number of particles in each direction.
     * @param m The mass of each particle.
     * @param distance The distance between neighbouring particles.
     * @param avgVelocityBrownian The average velocity of the Brownian motion.
     * @param dimension The dimension of the Brownian motion.
     * @param type The type of the particles.
     * @param epsilon The parameter epsilon of the Lennard-Jones potential.
     * @param sigma The parameter sigma of the Lennard-Jones potential.
     * @param stationary The flag that determines whether the cuboid is a wall.
     */
    Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
           double m, double distance, double avgVelocityBrownian, int dimension, int type, double epsilon, double sigma,
           bool stationary = false);

    /**
     * Create particles that belong to the cuboid and insert them into the given vector.
     * @param particles: The vector that particles of the cuboid should be added to.
     */
    void createParticles(std::vector<Particle>& particles) const override;
};


