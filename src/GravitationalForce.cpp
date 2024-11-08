#include <cmath>
#include "GravitationalForce.h"
#include "Particle.h"
#include "utils/vectorOperations.h"

GravitationalForce::GravitationalForce(double g): g(g) {}

GravitationalForce::GravitationalForce(): g(1) {}

std::array<double, 3> GravitationalForce::force(Particle particle1, Particle particle2) {
    return scalar_mul(
            this->g * particle1.getM() * particle2.getM() / pow(norm(sub(particle1.getX(), particle2.getX())), 3),
            sub(particle1.getX(), particle2.getX())
            );
}