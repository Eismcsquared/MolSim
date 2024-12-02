#include <array>
#include "body/Particle.h"
#include "container/BoundaryCondition.h"

#pragma once

class Force {
public:
    virtual std::array<double, 3> force(Particle& particle1, Particle& particle2) = 0;

    virtual ~Force() = default;
};
