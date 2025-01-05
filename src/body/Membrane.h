#include "Cuboid.h"

class Membrane : public Cuboid {

private:

    /**
     * The stiffness constant.
     */
    double k;

    /**
     * The average bond length.
     */
    double r0;

public:

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
     * @param k The stiffness constant of the bonding.
     * @param r0 The average bond length.
     */
    Membrane(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 2> &n, double m,
             double distance, double avgVelocityBrownian, int dimension, int type, double epsilon, double sigma,
             double k, double r0);
    /**
     * Create particles that belong to the membrane and insert them into the given vector.
     * @param particles: The vector that particles of the membrane should be added to.
     */
    void createParticles(std::vector<Particle>& particles) const override;

};
