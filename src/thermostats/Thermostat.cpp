#include <array>
#include "Thermostat.h"
#include "utils/ArrayUtils.h"



Thermostat::Thermostat(double targetT, int periode, double maxDelta, int dimension) : target_T(targetT), periode(periode),
                                                                                    maxDelta(maxDelta),
                                                                                    dimension(dimension) {}

void Thermostat::apply(std::vector<Particle> &particles) const {

    // compute the current temperature

    double T = temperature(particles, dimension);

    // If the current temperature is 0, the thermostat can't be applied.

    if (std::abs(T) < 1e-12) {
        return;
    }

    // compute the temperature after the application

    double TNew = T + 2 * ((target_T > T) - 0.5) * std::min(std::abs(target_T - T), maxDelta);

    // apply thermostat

    std::array<double, 3> meanVel = meanVelocity(particles);
    for (Particle &p: particles) {
        if (p.isInDomain()) {
            std::array<double, 3> newV = meanVel + sqrt(TNew / T) * (p.getV() - meanVel);
            p.setV(newV);
        }
    }
}

int Thermostat::getPeriode() const {
    return periode;
}

double Thermostat::getTargetT() const {
    return target_T;
}

double Thermostat::getMaxDelta() const {
    return maxDelta;
}

int Thermostat::getDimension() const {
    return dimension;
}


double Thermostat::temperature(std::vector<Particle> &particles, int dimension) {

    // compute the global mean velocity

    std::array<double, 3> meanVel = meanVelocity(particles);

    // compute the current temperature

    double E_kin = 0;
    int particleNumber = 0;
    for (Particle &p: particles) {
        if (p.isInDomain()) {
            E_kin += pow(ArrayUtils::L2Norm(p.getV() - meanVel), 2) / (2 * p.getM());
            particleNumber++;
        }
    }
    if (particleNumber == 0) {
        return 0;
    }
    return E_kin * 2 / (particleNumber * dimension);
}

std::array<double, 3> Thermostat::meanVelocity(std::vector<Particle> &particles) {
    std::array<double, 3> meanVelocity = {0, 0, 0};
    int particleNumber = 0;
    for (Particle &p: particles) {
        if (p.isInDomain()) {
            meanVelocity = meanVelocity + p.getV();
            particleNumber++;
        }
    }
    if (particleNumber == 0) {
        return {0, 0, 0};
    }

    return (1.0 / particleNumber) * meanVelocity;
}
