#include "Particle.h"

/**
 * @brief This class represents particles that do not move (e.g. for the wall).
 */
class StationaryParticle : public Particle{

    /**
     * Do nothing, as the particle is stationary.
     * @param deltaT The time step of the update.
     */
    virtual void updateX(double deltaT) override;

    /**
     * Do nothing, as the particle is stationary.
     * @param deltaT The time step of the update.
     */
    virtual void updateV(double deltaT) override;

    /**
     * Do nothing, as the particle is stationary.
     * @param force The force to be added as an array.
     */
    virtual void addForce(std::array<double, 3> force) override;
};

