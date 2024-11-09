#include <cmath>
#include "GravitationalForce.h"
#include "Particle.h"
#include "utils/vectorOperations.h"

GravitationalForce::GravitationalForce(double g): g(g) {}

GravitationalForce::GravitationalForce(): g(1) {}

std::array<double, 3> GravitationalForce::force(Particle particle1, Particle particle2) {
    std::array<double, 3> result{};
    double r = distance(particle1.getX(), particle2.getX());
    double factor = g * particle1.getM() * particle2.getM() / pow(r, 3);
    for (int i = 0; i < 3; ++i) {
        result[i] = factor * (particle1.getX()[i] - particle2.getX()[i]);
    }
    return result;
}

