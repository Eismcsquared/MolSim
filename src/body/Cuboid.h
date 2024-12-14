#include <array>
#include <vector>
#include "Particle.h"
#include "utils/ArrayUtils.h"
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
     * @param epsilon The parameter epsilon of the Lennard-Jones potential.
     * @param sigma The parameter sigma of the Lennard-Jones potential.
     */
    Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
           double m, double distance, double avgVelocityBrownian, int dimension, double epsilon, double sigma);

    /**
     * Create particles that belong to the cuboid and insert them into the given vector.
     * @param particles: The vector that particles of the cuboid should be added to.
     */
    inline void createParticles(std::vector<Particle>& particles) const override {
        for (unsigned int i = 0; i < N[0]; ++i) {
            for (unsigned int j = 0; j < N[1]; ++j) {
                for (unsigned int k = 0; k < N[2]; ++k) {
                    std::array<double, 3> x_particle = {x[0] + i * distance, x[1] + j * distance, x[2] + k * distance};
                    std::array<double, 3> vBrownian = maxwellBoltzmannDistributedVelocity(avgVelocityBrownian, dimension);
                    std::array<double, 3> v_particle = v + vBrownian;
                    particles.emplace_back(x_particle, v_particle, m, epsilon, sigma);
                }
            }
        }
    }
};


