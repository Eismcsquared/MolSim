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
    * @param start_time: The start time of the simulation. Default: 0.
    * @param end_time: The end time of the simulation.
    * @param delta_t: The time step of the simulation.
    * @param f: The force object that defines the force between two particles.
    * @param outputFormat: The output format of data, either .vtu or .xyz.
    */
   ParticleContainer(std::vector<Particle>& particles, double start_time, double end_time, double delta_t, Force* f,
    std::string outputFormat);


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
     * @brief Calculate the force between all particles version1
     */
    void updateF(bool newton3 = true);


     /**
      * @brief Update the position for all particles
      */
    void updateX();

    /**
     * @brief Update the velocity for all particles
     */
    void updateV();

    /**
     * @brief Calculate the position, force and velocity for all particles
     * @param version: 1 = without optimization or 2 = optimized for cache efficiency
     */
    void simulate(const std::string& out_name);

    /**
    * @brief Save the current state of the particles in the container to the output.
    * @param iteration: The current iteration number.
    */
    void plotParticles(int iterations, const std::string& out_name);

  
  private:


    /**
     * Particles stored in the container
     */
    std::vector<Particle>& particles;
    /**
     * Start time of the simulation
     */
    double start_time = 0;
    /**
     * End time of the simulation
     */
    double end_time = 1000;
    /**
     * Time step of the simulation
     */
    double delta_t;
    /**
     * The functional interface which computes the force between two particles (gravitational force, Lennard Jones force...)
     */
    Force& f;
    /**
     * Output format of the data, either .vtu or .xyz
     */
    std::string outputFormat;



};
