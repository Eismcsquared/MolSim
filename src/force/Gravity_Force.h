#include <cmath>
#include "body/Particle.h"
#include "utils/ArrayUtils.h"
#include "force/Force.h"


class Gravity_Force : public Force {
private:
    /**
     * gravitational acceleration
    */
    std::array<double, 3> g;

public:
    explicit Gravity_Force(const std::array<double, 3>& g);
    Gravity_Force() : g({0.0, -12.44, 0.0}) {}
    /**
     * gravitational force on gravitational field
     * particle2 is not used.
     */
    std::array<double, 3> force(Particle& particle1, Particle& particle2) override;

    ~Gravity_Force() override;
};