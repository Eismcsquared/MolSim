#include <array>
#include "utils/ArrayUtils.h"
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
    /**
     * Constructor. The parameters epsilon and sigma are set to 5 and 1 by default, respectively.
     * @param x The center of the sphere.
     * @param v The mean velocity of the sphere.
     * @param radius The number of particles along the radius.
     * @param m The mass of each particle.
     * @param distance The distance between neighbouring particles.
     * @param avgVelocityBrownian The average velocity of the Brownian motion.
     * @param dimension The dimension of the sphere.
     */
    Sphere(const std::array<double, 3> &x, const std::array<double, 3> &v, int radius,
           double m, double distance, double avgVelocityBrownian, int dimension);

    /**
     * Constructor. The parameters epsilon and sigma are set to 5 and 1 by default, respectively.
     * @param x The center of the sphere.
     * @param v The mean velocity of the sphere.
     * @param radius The number of particles along the radius.
     * @param m The mass of each particle.
     * @param distance The distance between neighbouring particles.
     * @param avgVelocityBrownian The average velocity of the Brownian motion.
     * @param dimension The dimension of the sphere.
     * @param epsilon The parameter epsilon of the Lennard-Jones potential.
     * @param sigma The parameter sigma of the Lennard-Jones potential.
     */
    Sphere(const std::array<double, 3> &x, const std::array<double, 3> &v, int radius,
           double m, double distance, double avgVelocityBrownian, int dimension, double epsilon, double sigma);

    /**
     * Create particles that belong to the sphere and insert them into the given vector.
     * @param particles: The vector that particles of the cuboid should be added to.
     */
    inline void createParticles(std::vector<Particle>& particles) const override {
        for (int i = -radius; i <= radius; ++i) {
            for (int j = -radius; j <= radius; ++j) {
                if (dimension == 2) {
                    std::array<double, 3> pos = x + std::array<double, 3> {i * distance, j * distance, 0};
                    if(ArrayUtils::L2Norm(x - pos) <= radius * distance) {
                        particles.emplace_back(pos, v + maxwellBoltzmannDistributedVelocity(avgVelocityBrownian, dimension), m, epsilon, sigma);
                    }
                } else {
                    for (int k = -radius; k <= radius; ++k) {
                        std::array<double, 3> pos = x + std::array<double, 3> {i * distance, j * distance, k * distance};
                        if(ArrayUtils::L2Norm(x - pos) <= radius * distance) {
                            particles.emplace_back(pos, v + maxwellBoltzmannDistributedVelocity(avgVelocityBrownian, dimension), m, epsilon, sigma);
                        }
                    }
                }
            }
        }
    }
};
