#include <string>
#include <memory>
#include "body/Particle.h"

/**
 * @brief Define an interface for iteration over particles in a container.
 */
class Iterator {
public:
    /**
     * Destructor.
     */
    virtual ~Iterator() = default;
     /**
      * Determine the next particle in the container and update the iterator.
      * @return The next particle in the container
      */
    virtual Particle& next() = 0;
     /**
      * Determine whether there are further particles in the container.
      * @return True if the end of the container is not yet reached.
      */
    virtual bool hasNext() = 0;
};

/**
 * @brief The abstract class that represents a container for particles. Subclasses implement the concrete way to manage particles.
 */
class Container {
public:
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
     * Update the force of particles based on their positions.
     */
    virtual void updateF() = 0;
    /**
     * Determine the number of particles in a container.
     * @return The number of particles contained in the container.
     */
    virtual unsigned long getParticleNumber() const = 0;

    /**
     * Return a new iterator for the container.
     * @return a new iterator for the container.
     */
    virtual std::unique_ptr<Iterator> iterator() const = 0;
};



