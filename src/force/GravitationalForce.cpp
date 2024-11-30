#include "GravitationalForce.h" 


GravitationalForce::GravitationalForce(double g): g(g) {}

GravitationalForce::GravitationalForce(): g(1) {}

std::array<double, 3> GravitationalForce::force(Particle& particle1, Particle& particle2) {
    double r = ArrayUtils::L2Norm(particle1.getX() - particle2.getX());
    double factor = g * particle1.getM() * particle2.getM() / pow(r, 3);
    return factor * (particle1.getX() - particle2.getX());
}

std::array<double, 3> GravitationalForce::ghostforce(Particle& particle, std::array<double, 3> boundary, double distance, std::array<bool, 6> boundarycondition) {
    return {0,0,0}; // no ghost force for gravitational force
}

GravitationalForce::~GravitationalForce() = default;

