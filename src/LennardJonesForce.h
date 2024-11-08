
#include "Force.h"

class LennardJonesForce: public Force {
private:
    double epsilon;
    double sigma;
public:
    LennardJonesForce(double epsilon, double sigma);
    LennardJonesForce();

    /**
     * @brief calculate the Lennard-Jones force between two particles.
     * @param particle1 the first particle.
     * @param particle2 the second particle.
     * @return the Lennard-Jones force of the first particle acting on the second particle
     */
    std::array<double, 3> force(Particle particle1, Particle particle2) override;
};
