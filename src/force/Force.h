#include <array>
#include "body/Particle.h"
#include "container/BoundaryCondition.h"

#pragma once

/**
 * @brief This is an abstract class that define the function for force calculation between two objects.
 */
class Force {
public:
    virtual std::array<double, 3> force(Particle& particle1, Particle& particle2) = 0;

    virtual ~Force() = default;
};
