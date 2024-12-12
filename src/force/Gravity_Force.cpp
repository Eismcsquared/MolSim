#include "Gravity_Force.h"

Gravity_Force::Gravity_Force(const std::array<double, 3>& g) : g(g) {}

std::array<double, 3> Gravity_Force::force(Particle& particle1, Particle& particle2) {
    return particle1.getM() * g;
}

Gravity_Force::~Gravity_Force() = default;