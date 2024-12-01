#include "LennardJonesForce.h"
#include "spdlog/spdlog.h"

LennardJonesForce::LennardJonesForce(double epsilon, double sigma): epsilon(epsilon), sigma(sigma) {}
LennardJonesForce::LennardJonesForce(): epsilon(5), sigma(1) {}

std::array<double, 3> LennardJonesForce::force(Particle& particle1, Particle& particle2) {
    double r = ArrayUtils::L2Norm(particle1.getX() - particle2.getX());
    double a = pow(sigma / r, 6);
    double factor = 24 * epsilon / pow(r, 2) * (a - 2 * pow(a, 2));
    return factor * (particle1.getX() - particle2.getX());
}

std::array<double, 3> LennardJonesForce::ghostForce(Particle &particle, std::array<double, 3> boundary,
                                                    std::array<BoundaryCondition, 6> boundaryCondition) {
    std::array<double, 3> gForce = {0, 0, 0};
    std::array<double, 3> pos = particle.getX();
    std::array<double, 3> vel = particle.getV();
    double m = particle.getM();
    std::array<Particle, 6> ghostParticles;

    ghostParticles[0] = Particle(std::array<double, 3>{-pos[0], pos[1], pos[2]}, std::array<double, 3>{-vel[0], vel[1], vel[2]}, m);
    ghostParticles[1] = Particle(std::array<double, 3>{2 * boundary[0] - pos[0], pos[1], pos[2]}, std::array<double, 3>{-vel[0], vel[1], vel[2]}, m);
    ghostParticles[2] = Particle(std::array<double, 3>{pos[0], -pos[1], pos[2]}, std::array<double, 3>{vel[0], -vel[1], vel[2]}, m);
    ghostParticles[3] = Particle(std::array<double, 3>{pos[0], 2 * boundary[1] - pos[1], pos[2]}, std::array<double, 3>{vel[0], -vel[1], vel[2]}, m);
    ghostParticles[4] = Particle(std::array<double, 3>{pos[0], pos[1], -pos[2]}, std::array<double, 3>{vel[0], vel[1], -vel[2]}, m);
    ghostParticles[5] = Particle(std::array<double, 3>{pos[0], pos[1], 2 * boundary[2] - pos[2]}, std::array<double, 3>{vel[0], vel[1], -vel[2]}, m);

    for (int i = 0; i < 6; ++i) {
        if (boundaryCondition[i] == REFLECTING && ArrayUtils::L2Norm(pos - ghostParticles[i].getX()) <= sigma * pow(2, 1.0 / 6)) {
            gForce = gForce + force(ghostParticles[i], particle);
        }
    }

    return gForce;
}


LennardJonesForce::~LennardJonesForce() = default;