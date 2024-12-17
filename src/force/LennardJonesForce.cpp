#include "LennardJonesForce.h"
#include "spdlog/spdlog.h"


std::array<double, 3> LennardJonesForce::force(Particle& particle1, Particle& particle2) {
    double epsilon;
    if (std::abs(particle1.getEpsilon() - particle2.getEpsilon()) < 1e-12) {
        epsilon = particle1.getEpsilon();
    } else {
        epsilon = sqrt(particle1.getEpsilon() * particle2.getEpsilon());
    }
    double sigma = (particle1.getSigma() + particle2.getSigma()) / 2;

    std::array<double, 3> r_21 = particle1.getX() - particle2.getX();
    double rSquared = ArrayUtils::L2NormSquare(r_21);
    double a = pow(sigma * sigma / rSquared, 3);
    double factor = 24 * epsilon / rSquared * (a - 2 * a * a);
    return factor * r_21;
}


LennardJonesForce::~LennardJonesForce() = default;