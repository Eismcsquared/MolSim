
#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include "Particle.h"
#include <vector>
#include <array>

#pragma once


/**
 * @brief Class that manages a system of particles. It stores particles in a vector, which allows efficient random access. Moreover, the positions of the particles are (redundantly) stored in an extra vector. This should increase cache efficiency when calculating gravitational forces. This class provides functionality to update parameters of the particles using Velocity-Störmer-Verlet.
 */
class ParticleContainer {

  public:

   ParticleContainer(std::vector<Particle> particles, double start_time, double end_time, double delta_t,
    std::string outputformat);

  /**
   * @brief Destructor
   */
   ~ParticleContainer();

  
    /**
     * @brief Add a particle to the container
     * @param particle: The particle to add to the container
     */
    void addParticle(Particle particle);

    /**
     * @brief Calculate the force between all particles version1
     */
    void calculateF_v1();

    /**
     * @brief Calculate the force between all particles version2 , optimized version by memory access
     */
    void calculateF_v2();

     /**
      * @brief Update the position for all particles
      */
    void calculateX();

    /**
     * @brief Update the velocity for all particles
     */
    void calculateV();

    /**
     * @brief Calculate the position, force and velocity for all particles
     * @param version: 1 = without optimization or 2 = optimized for cache efficiency
     */
    void calculate(int version);

    /**
    * @brief Save the current state of the particles in the container to the output.
    * @param iteration: The current iteration number.
    */
    void plotParticles(int iterations);

  
  private:


    /**
     * Particles stored in the container
     */
    std::vector<Particle> particles;
    /**
     * Positions of the particles stored separately for more efficient access
     */
    std::vector<std::array<double,3>> positions;
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
     * Output format of the data, either .vtu or .xyz
     */
    std::string outputformat;


};

#endif //PARTICLECONTAINER_H
