#include "Force.h"
#include "utils/ArrayUtils.h"

/**
 * @brief This class computes interactions between particles within a membrane, which is harmonic for neighbouring particles and repulsive Lennard-Jones otherwise.
 */
class MembraneForce: public Force {

    /**
     * @brief compute the membrane force between two particles. This force is harmonic if both particles are neighbours
     * and Lennard-Jones otherwise. The force between non-neighboring particles is only applied if they are closed to each
     * other and the corresponding Lennard-Jones force becomes repulsive.
     * @param particle1 the first particle.
     * @param particle2 the second particle.
     * @return the membrane force of the first particle acting on the second particle
     */
    std::array<double, 3> force(Particle& particle1, Particle& particle2) override;

};

