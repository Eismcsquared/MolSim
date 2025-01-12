#include "Sphere.h"


Sphere::Sphere(const std::array<double, 3> &x, const std::array<double, 3> &v, int radius, double m,
               double distance, double avgVelocityBrownian, int dimension) :
               Cluster(x, v, m, distance, avgVelocityBrownian, dimension), radius(radius){}


Sphere::Sphere(const std::array<double, 3> &x, const std::array<double, 3> &v, int radius, double m, double distance,
               double avgVelocityBrownian, int dimension, int type, double epsilon, double sigma):
        Cluster(x, v, m, distance, avgVelocityBrownian, dimension, type, epsilon, sigma), radius(radius){}

void Sphere::createParticles(std::vector<Particle> &particles) const {
    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {
            int r = dimension == 2 ? 0 : radius;
            for (int k = -r; k <= r; ++k) {
                std::array<double, 3> pos = x + std::array<double, 3> {i * distance, j * distance, k * distance};
                if(ArrayUtils::L2NormSquare(x - pos) <= radius * radius * distance * distance) {
                    particles.emplace_back(pos, v + maxwellBoltzmannDistributedVelocity(avgVelocityBrownian, dimension), m, type, epsilon, sigma);
                }
            }
        }
    }
}


