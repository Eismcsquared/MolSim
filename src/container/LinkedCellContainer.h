#include <vector>
#include "container/Container.h"
#include "container/Cell.h"

/**
 * @brief This class represents a particle container that implements the linked cell algorithm.
 */
class LinkedCellContainer: public Container {
private:

    /**
     * The cells that the domain is divided into.
     */
    std::vector<Cell> cells = std::vector<Cell>();

public:

    /**
    * Construct a linked cell container.
    * @param particles: The particles to store.
    * @param f: The force object that defines the force between two particles.
    */
    LinkedCellContainer(std::vector<Particle> &particles, std::unique_ptr<Force> &f_ptr);
    /**
     * @brief Update the force between all particles.
     * @param newton3 The Newton's third law is applied if this flag is set.
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

