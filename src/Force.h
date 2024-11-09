#include <array>
#include "Particle.h"

class Force {
public:
    virtual std::array<double, 3> force(Particle particle1, Particle particle2) = 0;
};
