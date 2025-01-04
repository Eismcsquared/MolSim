#include "MembraneForce.h"
#include "LennardJonesForce.h"

std::array<double, 3> MembraneForce::force(Particle &particle1, Particle &particle2) {
    std::array<double, 3> r21 = particle1.getX() - particle1.getX();
    std::array<double, 3> normalizedR21 = (1 / ArrayUtils::L2Norm(r21)) * r21;
    const double sqrt2 = 1.414213562373095;
    const double cbrt2 = 1.259921049894873;
    // neighbouring particles are assumed to have the same stiffness constant and average bond length.
    double k = particle1.getK();
    double r0 = particle1.getR0();
    double sigma = (particle1.getSigma() + particle2.getSigma()) / 2;
    if (particle1.isNeighbour(particle2)) {
        return k * (r21 - r0 * normalizedR21);
    } else if (particle1.isDiagonalNeighbour(particle2)) {
        return k * (r21 - sqrt2 * r0 * normalizedR21);
    } else if (ArrayUtils::L2NormSquare(r21) <= cbrt2 * sigma * sigma) {
        LennardJonesForce lennardJonesForce;
        return lennardJonesForce.force(particle1, particle2);
    } else {
        return {0, 0, 0};
    }
}
