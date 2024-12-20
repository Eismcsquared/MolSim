#include "GravitationalForce.h" 


GravitationalForce::GravitationalForce(double g): g(g) {}

GravitationalForce::GravitationalForce(): g(1) {}

std::array<double, 3> GravitationalForce::force(Particle& particle1, Particle& particle2) {
    double r = ArrayUtils::L2Norm(particle1.getX() - particle2.getX());
    double factor = g * particle1.getM() * particle2.getM() / pow(r, 3);
    return factor * (particle1.getX() - particle2.getX());
}
GravitationalForce::~GravitationalForce() = default;

