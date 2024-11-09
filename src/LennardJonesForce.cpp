#include "LennardJonesForce.h"
#include "utils/vectorOperations.h"

LennardJonesForce::LennardJonesForce(double epsilon, double sigma): epsilon(epsilon), sigma(sigma) {}
LennardJonesForce::LennardJonesForce(): epsilon(5), sigma(1) {}

std::array<double, 3> LennardJonesForce::force(Particle particle1, Particle particle2) {
    std::array<double, 3> result{};
    double r = distance(particle1.getX(), particle2.getX());
    double a = pow(sigma / r, 6);
    double factor = 24 * epsilon / pow(r, 2) * (a - 2 * pow(a, 2));
    for (int i = 0; i < 3; ++i) {
        result[i] = factor * (particle1.getX()[i] - particle2.getX()[i]);
    }
    return result;
}
