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

std::array<double, 3> LennardJonesForce::ghostforce(Particle& particle, std::array<double, 3> boundary, double distance, std::array<bool, 6> boundarycondition) {
    double h = distance / 2;
    std::array<double, 3> gforce = {0, 0, 0};
    std::array<double, 3> pos = particle.getX();

    auto calculateForce = [&](double coord, double boundaryLimit, bool condition, bool isLower) -> double {
        if (!condition) return 0.0;
        double r = std::abs(isLower ? 2 * coord : 2 * (boundaryLimit - coord));

        if (r >= h || r <= std::numeric_limits<double>::epsilon()) return 0.0;

        double a = pow(sigma / r, 6);
        double factor = 24 * epsilon / pow(r,2) * (a - 2 * pow(a, 2));
        return factor * (isLower ? 2 * coord : 2 * (coord - boundaryLimit));
    };

    for (int i = 0; i < gforce.size(); ++i) {
        if (boundary[i] != 0) { 
            gforce[i] += calculateForce(pos[i], boundary[i], boundarycondition[i * 2], true);  // Lower boundary
            gforce[i] += calculateForce(pos[i], boundary[i], boundarycondition[i * 2 + 1], false); // Upper boundary
        }
    }

    return gforce;
}


LennardJonesForce::~LennardJonesForce() = default;