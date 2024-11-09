
#include "Cuboid.h"
#include "utils/MaxwellBoltzmannDistribution.h"

std::vector<Particle> Cuboid::createParticles() {
    std::vector<Particle> particles;
    for (int i = 0; i < N[0]; ++i) {
        for (int j = 0; j < N[1]; ++j) {
            for (int k = 0; k < N[2]; ++k) {
                std::array<double, 3> x_particle = {x[0] + i * distance, x[1] + j * distance, x[2] + k * distance};
                std::array<double, 3> vBrownian = maxwellBoltzmannDistributedVelocity(avgVelocityBrownian, 3);
                std::array<double, 3> v_particle = {v[0] + vBrownian[0], v[1] + vBrownian[1], v[2] + vBrownian[2]};
                particles.emplace_back(x_particle, v_particle, m);
            }
        }
    }
    return particles;
}