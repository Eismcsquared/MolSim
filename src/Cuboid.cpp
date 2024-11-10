#include "utils/ArrayUtils.h"
#include "Cuboid.h"
#include "utils/MaxwellBoltzmannDistribution.h"

std::vector<Particle> Cuboid::createParticles() const {
    std::vector<Particle> particles = std::vector<Particle>();
    for (unsigned int i = 0; i < N[0]; ++i) {
        for (unsigned int j = 0; j < N[1]; ++j) {
            for (unsigned int k = 0; k < N[2]; ++k) {
                std::array<double, 3> x_particle = {x[0] + i * distance, x[1] + j * distance, x[2] + k * distance};
                std::array<double, 3> vBrownian = maxwellBoltzmannDistributedVelocity(avgVelocityBrownian, 3);
                std::array<double, 3> v_particle = v + vBrownian;
                particles.emplace_back(x_particle, v_particle, m);
            }
        }
    }
    return particles;
}

Cuboid::Cuboid(const std::array<double, 3> &x, const std::array<double, 3> &v, const std::array<unsigned int, 3> &n,
               double m, double distance, double avgVelocityBrownian) : x(x), v(v), N(n), m(m), distance(distance),
                                                                        avgVelocityBrownian(avgVelocityBrownian) {}
