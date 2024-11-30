#include <cmath>
#include "body/Particle.h"
#include "utils/ArrayUtils.h"
#include "force/Force.h"
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


    /**
     * @brief simulate the gravitational force between a particle and a ghost particle. 
     * @param particle the particle.
     * @param boundary the boundary of the system.
     * @param distance lennard jones potential distance
     * @param boundarycondition the boundary condition of the system.
     */
    std::array<double, 3> ghostforce(Particle& particle,std::array<double, 3> boundary, double distance,
    std::array<bool,6> boundarycondition) override;

    ~GravitationalForce() override;
};



