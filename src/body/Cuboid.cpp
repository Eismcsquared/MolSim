#include "Cuboid.h"

Cuboid::Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
               double m, double distance, double avgVelocityBrownian, int dimension) :
               Cluster(x, v, m, distance, avgVelocityBrownian, dimension), N(n) {}

Cuboid::Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
               double m, double distance, double avgVelocityBrownian, int dimension, int type, double epsilon, double sigma):
        Cluster(x, v, m, distance, avgVelocityBrownian, dimension, type, epsilon, sigma), N(n){}

void Cuboid::createParticles(std::vector<Particle> &particles) const {
    for (unsigned int i = 0; i < N[0]; ++i) {
        for (unsigned int j = 0; j < N[1]; ++j) {
            for (unsigned int k = 0; k < N[2]; ++k) {
                std::array<double, 3> x_particle = {x[0] + i * distance, x[1] + j * distance, x[2] + k * distance};
                std::array<double, 3> vBrownian = maxwellBoltzmannDistributedVelocity(avgVelocityBrownian, dimension);
                std::array<double, 3> v_particle = v + vBrownian;
                particles.emplace_back(x_particle, v_particle, m, type, epsilon, sigma);
            }
        }
    }
}


