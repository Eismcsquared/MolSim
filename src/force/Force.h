#include <array>
#include "body/Particle.h"

#pragma once

class Force {
public:
    virtual std::array<double, 3> force(Particle& particle1, Particle& particle2) = 0;

    virtual std::array<double, 3> ghostforce(Particle& particle,std::array<double, 3> boundary, double distance,
    std::array<bool,6> boundarycondition) = 0;

    virtual ~Force() = default;
};
