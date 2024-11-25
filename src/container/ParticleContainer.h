#include <vector>
#include <array>
#include <memory>
#include "body/Particle.h"
#include "body/Cuboid.h"
#include "force/Force.h"
#include "container/Container.h"

#pragma once

/**
 * @brief The Iterator for iteration over particles in a particle container.
 */
class ParticleContainerIterator: public Iterator {
    /**
     * The current particle.
     */
    std::vector<Particle>::iterator current;
    /**
     * The end of the iteration.
     */
    std::vector<Particle>::iterator end;
public:
    /**
     * Destructor.
     */
    ~ParticleContainerIterator() override = default;
    /**
     * Constructor.
     * @param current
     * @param end
     */
    ParticleContainerIterator(std::vector<Particle>::iterator current, std::vector<Particle>::iterator end);
    /**
     * Update the iterator and return the current particle.
     * @return The current particle.
     */
    Particle& next() override;
    /**
     * Determine whether there are further particles.
     * @return True if the end of iteration is not yet reached.
     */
    bool hasNext() override;
};

/**
 * @brief The class that implements the direct sum method. The parameters of the particles are updated using Velocity-Störmer-Verlet.
 */
class ParticleContainer: Container{

  public:

   /**
    * @brief Construct a particle container.
    * @param particles: The particles to store.
    * @param f: The force object that defines the force between two particles.
    */
   ParticleContainer(std::vector<Particle>& particles, std::unique_ptr<Force>& f, bool netwon3);


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
     * @brief Update the force between all particles.
     */
    void updateF() override;


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

    /**
     * Determine the number of particles in a container.
     * @return The number of particles contained in the container.
     */
    unsigned long getParticleNumber() const override;

    /**
     * The getter for the particles in a particle container.
     * @return The particles in the container as a vector.
     */
    std::vector<Particle>& getParticles() const;

    /**
     * Provide a string representation of a container.
     * @return A string representation of the container.
     */
    std::string toString();

    /**
     * Compare the current container with another container based on the positions, velocities and forces of particles in the container.
     * @param other The other particle container that should be compared with.
     * @return True if both particle containers contain the same particles.
     */
    bool operator==(const ParticleContainer& other) const;

    /**
     * Return a new iterator for the particle container.
     * @return A new iterator for the particle container.
     */
    std::unique_ptr<Iterator> iterator() const override;


private:


    /**
     * Particles stored in the container
     */
    std::vector<Particle>& particles;

    /**
     * The functional interface which computes the force between two particles (gravitational force, Lennard Jones force...)
     */
    Force& f;

    /**
     * Optional parameter to specify whether Newton's third law should be applied in the force calculation.
     */
    bool newton3;
};
