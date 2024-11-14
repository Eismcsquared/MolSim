#include "Particle.h"
#include <vector>
#include <array>
#include "Cuboid.h"
#include "Force.h"

#pragma once


/**
 * @brief Class that manages a system of particles. It stores particles in a vector, which allows efficient random access. Moreover, the positions of the particles are (redundantly) stored in an extra vector. This should increase cache efficiency when calculating gravitational forces. This class provides functionality to update parameters of the particles using Velocity-Störmer-Verlet.
 */
class ParticleContainer {

  public:

    /**
    * @brief Construct a particle container.
    * @param particles: The particles to store.
    * @param f: The force object that defines the force between two particles.
    */
   ParticleContainer(std::vector<Particle>& particles, Force* f);


/**
   * @brief Destructor
   */
   ~ParticleContainer();

  
    /**
     * @brief Add a particle to the container
     * @param particle: The particle to add to the container
     */
    void addParticle(const Particle& particle);

    /**
     * @brief Add the particles in a cuboid to the container
     * @param cuboid: The cuboid to add to the container
     */
    void addCuboid(const Cuboid &cuboid);

    /**
     * @brief Calculate the force between all particles
     * @param newton3: Specifies whether the Newton's third law should be applied in the force calculation.
     */
    void updateF(bool newton3 = true);


     /**
      * @brief Update the position for all particles
      * @param delta_t: The duration that the positions should be updated for.
      */
    void updateX(double delta_t);

    /**
     * @brief Update the velocity for all particles
     * @param delta_t: The duration that the velocities should be updated for.
     */
    void updateV(double delta_t);

    /**
     * @brief Simulate the system.
     * @param delta_t: The length of each time step.
     * @param end_time: The duration of the simulation.
     * @param out_name: The name of the files that data should be written to.
     * @param output_format: The format of the output, should be either vtu or xyz.
     */
    void simulate(double delta_t, double end_time, const std::string& out_name, const std::string& output_format, bool save_output = true);

    /**
     * @brief Save the current state of the particles in the container to the output.
     * @param iteration: The current iteration number.
     * @param out_name: The name of the files that data should be written to.
     * @param output_format: The format of the output, should be either vtu or xyz.
    */
    void plotParticles(int iterations, const std::string& out_name, const std::string& output_format);

    int getParticleNumber() const;

    std::vector<Particle>& getParticles() const;

    std::string toString();

    bool operator==(const ParticleContainer& other) const;

  
  private:


    /**
     * Particles stored in the container
     */
    std::vector<Particle>& particles;

    /**
     * The functional interface which computes the force between two particles (gravitational force, Lennard Jones force...)
     */
    Force& f;

};
