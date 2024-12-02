#include <vector>
#include "container/ParticleContainer.h"
#include "container/Cell.h"
#include "container/BoundaryCondition.h"

/**
 * @brief represents a direction in 3D
 */
enum Direction {
    LEFT, // negative x-direction
    RIGHT, // positive x-direction
    DOWN, // negative y-direction
    UP, // positive y-direction
    BACK, //negative z-direction
    FRONT // positive z-direction
};
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
     * The number of cells in each dimension within the domain.
     */
    std::array<int, 3> nCells;

    /**
     * The boundary condition of the domain, stored in the order: left, right, down, up, back, front
     * where: left->right is the x-direction, down->up is the y-direction, back->front is the z-direction.
     */
    std::array<BoundaryCondition, 6> boundaryConditions;

 

public:

    /**
    * Construct a linked cell container.
    * @param particles: The particles to store.
    * @param f: The force object that defines the force between two particles.
    */
    LinkedCellContainer(std::unique_ptr<std::vector<Particle>>& particles, std::unique_ptr<Force>& f,
                        std::array<double, 3> domainSize,double cutoff, std::array<BoundaryCondition, 6> boundaryConditions);
    /**
     *  @brief Update the force between all particles.
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
    int getCellIndex(std::array<double, 3> positions);

    /**
     * @brief Convert 3D cell index to 1D cell index.
     * @param The 3D cell index.
     * @return The 1D cell index.
     */
    int get1DIndex(std::array<int, 3> index3D);

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
    void updateCellF(const std::vector<int> &v1, const std::vector<int> &v2, bool newton3);

    /**
     * Check whether a cell is a boundary cell.
     * @param index The linear index of the cell.
     * @return True if the cell is at a boundary.
     */
    bool isBoundaryCell(int index);

    /**
     * Check whether a cell is a halo cell.
     * @param index The linear index of the cell.
     * @return True if the cell is a halo cell.
     */
    bool isHaloCell(int index);

    /**
     * Compute the indices of all halo cells
     * @param direction: The direction of the boundary, e.g. for LEFT: Indices of all halo cells with x < 0 are returned
     * @return The indices of halo cells in a vector
     */
    std::vector<int> getAllHaloIndices(Direction direction);

    /**
     * Remove particles from halo cells.
     * @param direction The direction in which the halo cells should be cleared.
     */
    void removeFormHalo(Direction direction);

    /**
     * Handle particles in the halo cells in one direction according to the boundary condition, i.e.
     * For OUTFLOW boundary condition: particles are removed form the domain.
     * For REFLECTING boundary condition: the corresponding component of the velocity is inverted.
     * @param direction The direction of the halo cells
     * @param boundaryCondition The boundary condition.
     */
    void updateHalo(Direction direction, BoundaryCondition boundaryCondition);

    /**
     * Handle particles in the halo cells according to the boundary condition.
     */
     void updateHalo();

    /**
     * Add a particle to the linked cell container
     * @param particle: The particle to add to the container
     */
    void addParticle(const Particle& particle) override;

    /**
     * Add the particles in a cluster to the container
     * @param cluster: The cluster to add to the container
     */
    void addCluster(const Cluster &cluster) override;

    bool operator==(const LinkedCellContainer &other) const;
};

