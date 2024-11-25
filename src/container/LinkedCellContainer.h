#include <vector>
#include "container/Container.h"
#include "container/Cell.h"

/**
 * @brief The iterator that iterates over the particles in a linked cell container.
 */
class LinkedCellContainerIterator: public Iterator {
    /**
     * The current particle.
     */
    std::vector<Particle>::iterator currentParticle;
    /**
     * The end of the current cell.
     */
    std::vector<Particle>::iterator endCurrentCell;

    /**
     * The current cell.
     */
    std::vector<Cell>::iterator currentCell;

    /**
     * The end of the iteration.
     */
    std::vector<Cell>::iterator end;
public:
    /**
     * Destructor.
     */
    ~LinkedCellContainerIterator() override = default;
    /**
     * Constructor.
     * @param current The firste cell.
     * @param end The end of cells.
     */
    LinkedCellContainerIterator(std::vector<Cell>::iterator currentCell, std::vector<Cell>::iterator end);
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
 * @brief This class represents a particle container that implements the linked cell algorithm.
 */
class LinkedCellContainer: Container {
private:

    /**
     * The cells that the domain is divided into.
     */
    std::vector<Cell>& cells;

public:
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
     * Determine the number of particles in a linked cell container.
     * @return The number of particles contained in the linked cell container.
     */
    unsigned long getParticleNumber() const override;

    /**
     * Return a new iterator for the linked cell container.
     * @return a new iterator for the linked cell container.
     */
    virtual std::unique_ptr<Iterator> iterator() const = 0;


};

