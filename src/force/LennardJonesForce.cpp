#include "LennardJonesForce.h"
#include "spdlog/spdlog.h"


std::array<double, 3> LennardJonesForce::force(Particle& particle1, Particle& particle2) {
    double epsilon = sqrt(particle1.getEpsilon() * particle2.getEpsilon());
    double sigma = (particle1.getSigma() + particle2.getSigma()) / 2;
    double r = ArrayUtils::L2Norm(particle1.getX() - particle2.getX());
    double a = pow(sigma / r, 6);
    double factor = 24 * epsilon / pow(r, 2) * (a - 2 * pow(a, 2));
    return factor * (particle1.getX() - particle2.getX());
}


LennardJonesForce::~LennardJonesForce() = default;