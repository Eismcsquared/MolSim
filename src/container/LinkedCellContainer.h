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
     * Store a pair of integers.
     */
    struct Pair{
        int first;
        int second;
    };

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
     * The indices of halo cells in each direction. Each vector represents a direction, following the order of the enum Direction.
     */
    std::array<std::vector<int>, 6> haloCells;

    /**
     * The indices of domain cells.
     */
    std::vector<int> domainCells;

    /**
     * The boundary condition of the domain, stored in the order: left, right, down, up, back, front
     * where: left->right is the x-direction, down->up is the y-direction, back->front is the z-direction.
     */
    std::array<BoundaryCondition, 6> boundaryConditions;

    /**
     * Store all pairs of neighbouring cells in vectors, in which every cell appears at most once.
     */
    std::vector<std::vector<Pair>> cellPairs;

    /**
     * Initialize the pairs of cells.
     */
    void initializePairs();
public:

     /**
     * Construct a linked cell container.
     * @param particles: The particles to store.
     * @param f: The force object that defines the force between two particles.
     * @param domainSize: The size of the domain.
     * @param cutoff: The cutoff radius.
     * @param boundaryConditions: The boundary conditions on different boundaries.
     */
     LinkedCellContainer(std::vector<Particle>& particles, std::unique_ptr<Force> &f,
                         std::array<double, 3> domainSize, double cutoff, std::array<BoundaryCondition, 6> boundaryConditions);


    /**
     * Getter for the cells.
     * @return The cells.
     */
    const std::vector<Cell> &getCells() const;

    /**
     * Getter for the cutoff radius.
     * @return The cutoff radius.
     */
    double getCutoff() const;

    /**
     * Getter for the domain size.
     * @return The domain size.
     */
    const std::array<double, 3> &getDomainSize() const;

    /**
     * Getter for the number of cells in the domain.
     * @return The number of cells in the domain.
     */
    const std::array<int, 3> &getNCells() const;

    /**
     * Getter for the boundary conditions.
     * @return The boundary conditions.
     */
    const std::array<BoundaryCondition, 6> &getBoundaryConditions() const;

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
     * @brief Get the neighbor cells of a cell.
     * @param cellIndex The index of the cell.
     * @return The indices of the neighbor cells.
     */
    std::set<int> getNeighborCells(int cellIndex);

    /**
     * @brief Get the index of a cell that contains a given position.
     * @param positions The position that should be tested.
     * @return The index of the cell that contains the given position.
     */
    inline int getCellIndex(std::array<double, 3> positions) {

        std::array<double, 3> cellSize = cells[0].getSize();
        // Check if the position is within the domain + halo region
        if(positions[0] < -cellSize[0] || positions[0] >= domainSize[0] + cellSize[0] || positions[1] < -cellSize[1] ||
           positions[1] >= domainSize[1] + cellSize[1] || positions[2] < -cellSize[2] || positions[2] >= domainSize[2] + cellSize[2]){
            return -1;
        }

        int idxX = std::floor(positions[0] / cellSize[0]);
        int idxY = std::floor(positions[1] / cellSize[1]);
        int idxZ = std::floor(positions[2] / cellSize[2]);

        // Calculate the linear cell index
        int cellIndex = get1DIndex({idxX, idxY, idxZ});

        return cellIndex;
    }

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
     * @param c1 The index of the first cell.
     * @param v2 The index of the second cell.
     */
    void updateFCells(int c1, int c2);

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
     * Check whether a cell is in domain.
     * @param index The linear index of the cell.
     * @return True if the cell is in domain.
     */
    bool isDomainCell(int index);

    /**
     * Remove particles from halo cells.
     * @param direction The direction in which the halo cells should be cleared.
     */
    void removeFromHalo(Direction direction);

    /**
     * Handle particles in the halo cells in one direction according to the boundary condition, i.e.
     * For OUTFLOW boundary condition: particles are removed form the domain.
     * For REFLECTING boundary condition: the corresponding component of the velocity is inverted.
     * @param direction The direction of the halo cells
     * @param boundaryCondition The boundary condition.
     * @param deltaT The time step.
     */
    void updateHalo(Direction direction, BoundaryCondition boundaryCondition, double deltaT);

    /**
     * Handle particles in the halo cells according to the boundary condition.
     */
     void updateHalo(double deltaT);

    /**
     * Add a particle to the linked cell container
     * @param particle: The particle to add to the container
     * @param deltaT The time step.
     */
    void addParticle(const Particle& particle) override;

    /**
     * Add the particles in a cluster to the container
     * @param cluster: The cluster to add to the container
     */
    void addCluster(const Cluster &cluster) override;

    /**
     * Provide a string representation of the linked cell container.
     * @return A string representation of the container.
     */
    std::string toString() override;

    bool operator==(const LinkedCellContainer &other) const;
};

