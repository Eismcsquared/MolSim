#include "utils/ArrayUtils.h"
#include "force/Force.h"

/**
 * @brief Class for calculation of the force between two particles assuming a Lennard-Jones potential.
 */
class LennardJonesForce: public Force {
public:

    /**
     * @brief simulate the Lennard-Jones force between two particles.
     * @param particle1 the first particle.
     * @param particle2 the second particle.
     * @return the Lennard-Jones force of the first particle acting on the second particle
     */
    std::array<double, 3> force(Particle& particle1, Particle& particle2) override;


    ~LennardJonesForce() override;
};

