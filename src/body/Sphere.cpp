#include "Sphere.h"
#include "utils/ArrayUtils.h"

Sphere::Sphere(const std::array<double, 3> &x, const std::array<double, 3> &v, int radius, double m,
               double distance, double avgVelocityBrownian, int dimension) :
               Cluster(x, v, m, distance, avgVelocityBrownian, dimension), radius(radius){}

void Sphere::createParticles(std::vector<Particle> &particles) const {
    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {
            if (dimension == 2) {
                std::array<double, 3> pos = x + std::array<double, 3> {i * distance, j * distance, 0};
                if(ArrayUtils::L2Norm(x - pos) <= radius * distance) {
                    particles.emplace_back(pos, v + maxwellBoltzmannDistributedVelocity(avgVelocityBrownian, dimension), m);
                }
            } else {
                for (int k = -radius; k <= radius; ++k) {
                    std::array<double, 3> pos = x + std::array<double, 3> {i * distance, j * distance, k * distance};
                    if(ArrayUtils::L2Norm(x - pos) <= radius * distance) {
                        particles.emplace_back(pos, v + maxwellBoltzmannDistributedVelocity(avgVelocityBrownian, dimension), m);
                    }
                }
            }
        }
    }
}