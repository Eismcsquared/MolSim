#include <string>
#include <memory>
#include <spdlog/spdlog.h>
#include <atomic>
#include "body/Particle.h"
#include "body/Cuboid.h"
#include "force/Force.h"
#include "thermostats/Thermostat.h"
#include "statistics/Statistics.h"

#pragma once

/**
 * @brief Define an interface for iteration over particles in a container.
 */
class Iterator {
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
     * Constructor.
     * @param begin The first particle.
     * @param end The end of the iteration.
     */
    Iterator(std::vector<Particle>::iterator begin, std::vector<Particle>::iterator end);
    /**
     * Destructor.
     */
    ~Iterator() = default;
     /**
      * Determine the next particle in the container and update the iterator.
      * @return The next particle in the container
      */
    Particle& next();
     /**
      * Determine whether there are further particles in the container.
      * @return True if the end of the container is not yet reached.
      */
    bool hasNext();
};

/**
 * @brief Represents external forces acting on single particles until a certain time.
 */
struct ConstantForce{
    unsigned int particleIndex;
    std::array<double, 3> f;
    double untilTime;
};

/**
 * @brief The abstract class that represents a container for particles. Subclasses implement the concrete way to manage particles.
 */
class ParticleContainer {

protected:
    /**
     * Particles stored in the container
     */
    std::vector<Particle> &particles;

    /**
     * The functional interface which computes the forces between two particles (gravitational force, Lennard Jones force...)
     */
    std::unique_ptr<Force> force;

    /**
     * The gravitational acceleration acting on all particles.
     */
    std::array<double, 3> g;

     /**
      * The number of particles in domain.
      */
    std::atomic<unsigned long> particleNumber;

     /**
      * The thermostat that adjusts the temperature of the system.
      */
    std::unique_ptr<Thermostat> thermostat;

    /**
     * The current time.
     */
    double t;

    /**
     * External forces acting on single particles until a certain time.
     */
    std::vector<ConstantForce> externalForces;

public:
    /**
    * Construct a particle container.
    * @param particles: The particles to store.
    * @param f_ptr: The force objects that defines the force between two particles.
    */
    ParticleContainer(std::vector<Particle> &particles, std::unique_ptr<Force> &f_ptr);

    /**
    * Construct a particle container.
    * @param particles: The particles to store.
    * @param f_ptr: The force objects that defines the force between two particles.
    * @param g: The gravitational acceleration.
    */
    ParticleContainer(std::vector<Particle> &particles, std::unique_ptr<Force> &f_ptr, std::array<double, 3> g);

    virtual ~ParticleContainer() = default;

    /**
     * The getter for the particles in a particle container that are still in domain.
     * @return The particles in the container as a vector.
     */
    std::vector<Particle>& getParticles() const;

    /**
     *
     * @return
     */
    Force &getForce() const;

    /**
     * Getter for the gravitational acceleration.
     * @return The value of g.
     */
    std::array<double, 3> getG() const;

    /**
     * Getter for the thermostat.
     * @return The thermostat.
     */
    std::unique_ptr<Thermostat> &getThermostat();

    /**
     * Determine the number of particles in a container.
     * @return The number of particles contained in the container.
     */
    unsigned long getParticleNumber() const;

    /**
     * The getter for the current time.
     * @return The current time.
     */
    double getT() const;

    /**
     * The getter for the external forces.
     * @return The external forces.
     */
    const std::vector<ConstantForce> &getExternalForces() const;
    /**
     * Setter for the force.
     * @param f
     */
    void setForce(std::unique_ptr<Force> &f);

    /**
     * Setter for the gravitational acceleration.
     * @param g
     */
    void setG(std::array<double, 3> g);

    /**
     * Setter for the thermostat.
     * @param thermostat A unique pointer to the thermostat.
     */
    void setThermostat(std::unique_ptr<Thermostat> &thermostat);

    /**
     * Setter for the time.
     * @param t The new time.
     */
    void setT(double t);

    /**
     * Update the position of particles by a time step.
     * @param delta_t The length of a time step.
     */
    virtual void updateX(double delta_t) = 0;
    /**
     * Update the velocity of particles by a time step.
     * @param delta_t The length of a time step.
     */
    virtual void updateV(double delta_t) = 0;


    /**
     * @brief Update the force between all particles using the default parallelization strategy.
     */
    inline void updateF() {updateF(0);}

    /**
     * Update the force between all particles.
     * @param strategy: The parallelization strategy.
     */
    virtual void updateF(int strategy) = 0;

    /**
     * Set the old force of particles to their current force and reset their force to the sum of external forces.
     */
    void resetF();

    /**
     * Add a particle to the container
     * @param particle: The particle to add to the container
     */
    virtual void addParticle(const Particle& particle);

    /**
     * Add the particles in a cluster to the container
     * @param cluster: The cluster to add to the container
     */
    virtual void addCluster(const Cluster &cluster);

    /**
     * Add an external force on a specific particle.
     * @param particleIndex The index of the particle the force should be applied on.
     * @param f The force as an array.
     * @param until The time until which the force is applied.
     */
    void addExternalForce(unsigned int particleIndex, const std::array<double, 3> f, double until);

    /**
     * Simulate the system of particles.
     * @param end_time The duration of the simulation.
     * @param delta_t The time step of the simulation.
     * @param out_name The name of the output file.
     * @param output_format The format of the output file, either "vtu" or "xyz".
     * @param output_frequency THe frequency of the output, in number of time steps.
     * @param save_output Output is activated if this flag is set.
     * @param strategy The parallelization strategy.
     * @param statistics Optional. If specified, statistics about the system is computed and stored.
     */
    void simulate(double end_time, double delta_t, const std::string& out_name,
                          const std::string& output_format, unsigned int output_frequency,
                          bool save_output, int strategy = 0, const std::shared_ptr<Statistics>& statistics = nullptr);

    /**
     * Write the current state of the container to the output files.
     * @param iteration The current iteration
     * @param out_name .The name of the output file.
     * @param output_format The format of the output file, either "vtu" or "xyz".
     */
    void plotParticles(int iteration, const std::string& out_name, const std::string& output_format);

    /**
     * Provide a string representation of a container.
     * @return A string representation of the container.
     */
    virtual std::string toString();

    /**
     * Compare the current container with another container based on the positions, velocities and forces of particles in the container.
     * @param other The other particle container that should be compared with.
     * @return True if both particle containers contain the same particles.
     */
    bool operator==(const ParticleContainer& other) const;

    /**
     * Return a new iterator for the container.
     * @return a new iterator for the container.
     */
    std::unique_ptr<Iterator> iterator() const;


};



