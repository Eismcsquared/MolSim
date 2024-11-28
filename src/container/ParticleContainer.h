#include <vector>
#include <array>
#include <memory>
#include "body/Particle.h"
#include "force/Force.h"
#include "container/Container.h"

#pragma once

/**
 * @brief The class that implements the direct sum method. The parameters of the particles are updated using Velocity-Störmer-Verlet.
 */
class ParticleContainer: public Container{

  public:

   /**
    * Construct a particle container.
    * @param particles: The particles to store.
    * @param f: The force object that defines the force between two particles.
    */
   ParticleContainer(std::unique_ptr<std::vector<Particle>>& particles, std::unique_ptr<Force>& f);


   /**
    * Destructor.
    */
   ~ParticleContainer();


    /**
     * Update the force between all particles.
     */
    void updateF(bool newton3) override;


    /**
     * @brief Update the position for all particles.
     * @param delta_t: The duration that the positions should be updated for.
     */
    void updateX(double delta_t) override;

    /**
     * @brief Update the velocity for all particles.
     * @param delta_t: The duration that the velocities should be updated for.
     */
    void updateV(double delta_t) override;


};
