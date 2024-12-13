#include <array>
#include "Thermostat.h"
#include "utils/ArrayUtils.h"



Thermostat::Thermostat(double targetT, int periode, double maxDelta, int dimension) : target_T(targetT), periode(periode),
                                                                                    maxDelta(maxDelta),
                                                                                    dimension(dimension) {}

void Thermostat::apply(std::vector<Particle> &particles) const {

    // compute the global mean velocity

    std::array<double, 3> meanVelocity = {0, 0, 0};
    int particleNumber = 0;
    for (Particle &p: particles) {
        if (p.isInDomain()) {
            meanVelocity = meanVelocity + p.getV();
            particleNumber++;
        }
    }
    meanVelocity = 1 / particleNumber * meanVelocity;

    // compute the current temperature

    double E_kin = 0;
    for (Particle &p: particles) {
        if (p.isInDomain()) {
            E_kin += pow(ArrayUtils::L2Norm(p.getV() - meanVelocity), 2) / (2 * p.getM());
        }
    }
    double temperature = E_kin * 2 / (particleNumber * dimension);

    // compute the temperature after the application

    double temperatureNew = temperature + 2 * ((target_T > temperature) - 0.5) * std::min(std::abs(target_T - temperature), maxDelta);

    // apply thermostat

    for (Particle &p: particles) {
        if (p.isInDomain()) {
            std::array<double, 3> newV = meanVelocity + sqrt(temperatureNew / temperature) * (p.getV() - meanVelocity);
            p.setV(newV);
        }
    }
}

int Thermostat::getPeriode() const {
    return periode;
}
