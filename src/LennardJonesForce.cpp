#include "LennardJonesForce.h"
#include "utils/vectorOperations.h"

LennardJonesForce::LennardJonesForce(double epsilon, double sigma): epsilon(epsilon), sigma(sigma) {}
LennardJonesForce::LennardJonesForce(): epsilon(5), sigma(1.1225) {}

std::array<double, 3> LennardJonesForce::force(Particle particle1, Particle particle2) {
    double distance = norm(sub(particle1.getX(), particle2.getX()));
    return {0, 0, 0};
}