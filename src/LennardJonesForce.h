
#include "Force.h"

/**
 * @brief Class for calculation of the force between two particles assuming a Lennard-Jones potential.
 */
class LennardJonesForce: public Force {
private:
    /**
     * @brief The parameter epsilon in the Lennard-Jones potential. Default: 5.
     */
    double epsilon;
    /**
     * @brief The parameter sigma in the Lennard-Jones potential. Default: 1.
     */
    double sigma;
public:
    LennardJonesForce(double epsilon, double sigma);
    LennardJonesForce();

    /**
     * @brief simulate the Lennard-Jones force between two particles.
     * @param particle1 the first particle.
     * @param particle2 the second particle.
     * @return the Lennard-Jones force of the first particle acting on the second particle
     */
    virtual std::array<double, 3> force(Particle& particle1, Particle& particle2) override;

    virtual ~LennardJonesForce() override;
};

