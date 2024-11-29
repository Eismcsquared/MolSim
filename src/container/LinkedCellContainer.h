#include <vector>
#include "container/ParticleContainer.h"
#include "container/Cell.h"

/**
 * @brief This class represents a particle container that implements the linked cell algorithm.
 */
class LinkedCellContainer: public ParticleContainer {
private:

    /**
     * The cells that the domain is divided into.
     */
    std::vector<Cell> cells = std::vector<Cell>();

    /**
     * The cutoff distance for the linked cell algorithm.
     */
    double cutoff;

    /**
     * The domain size of the simulation.
     */
    std::array<double, 3> domainSize;

    /**
     * The number of cells in each dimension.
     */
    std::array<int, 3> nCells;

    /**
     * The flag to determine if the domain should be reflected.
     */
    bool reflect = true;

    /**
     * The flag to determine if the particle should affected by ghost particles.
     */
    std::vector<bool> near_boundary;

    

public:

    /**
    * Construct a linked cell container.
    * @param particles: The particles to store.
    * @param f: The force object that defines the force between two particles.
    */
    LinkedCellContainer(std::unique_ptr<std::vector<Particle>>& particles, std::unique_ptr<Force>& f, std::array<double, 3> domainSize,double cutoff, bool reflect);
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
    
    /**
     * @brief Get the neighbor cells of a cell.
     * @param cellIndex The index of the cell.
     * @return The indices of the neighbor cells.
     */
    std::vector<int> getNeighborCells(int cellIndex);

    /**
     * @brief Get the index of a cell that contains a given position.
     * @param positions The position that should be tested.
     * @return The index of the cell that contains the given position.
     */
    unsigned int getCellIndex(std::array<double, 3> positions);

    /**
     * @brief Convert a 1D cell index to a 3D cell index.
     * @param cellIndex The 1D cell index.
     * @return The 3D cell index.
     */
    std::array<int, 3> get3DIndex(int cellIndex);


    /**
     * @brief Update the forces between particles in two cells.
     * @param v1 The indices of particles in the first cell.
     * @param v2 The indices of particles in the second cell.
     */
    void updateCellF(const std::vector<unsigned int> &v1, const std::vector<unsigned int> &v2, bool newton3);
};

