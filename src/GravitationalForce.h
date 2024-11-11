
#include "Force.h"
/**
 * @brief Class for calculation of the gravitational force between two objects.
 */
class GravitationalForce: public Force {
private:
    /**
     * the Newton constant, set to 1 by default
     */
    double g;
public:
    explicit GravitationalForce(double g);
    GravitationalForce();

    /**
     * @brief simulate the gravitational force between two particles.
     * @param particle1 the first particle.
     * @param particle2 the second particle.
     * @return the gravitational force of the first particle acting on the second particle
     */
    std::array<double, 3> force(Particle& particle1, Particle& particle2) override;

    ~GravitationalForce() override;
};



