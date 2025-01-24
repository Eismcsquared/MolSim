#include "container/LinkedCellContainer.h"
#include "spdlog/spdlog.h"
#include <vector>
#include <array>
#include "utils/ArrayUtils.h"
#include "cmath"
#include <omp.h>    


#define THREADS 4

// Will compare the computation time with and without nested parallelism (maybe in linux cluster faster because it has more cores and large cache)
#define InnerThreads 1

LinkedCellContainer::LinkedCellContainer(std::vector<Particle>& particles, std::unique_ptr<Force> &f,
                                         std::array<double, 3> domainSize, double cutoff, std::array<BoundaryCondition, 6> boundaryConditions) :
        ParticleContainer(particles, f), cutoff(cutoff), domainSize(domainSize), boundaryConditions(boundaryConditions){

    // set nested parallelism 
    omp_set_max_active_levels(1);

    // number of cells in each direction, making a cell larger if the cutoff radius does not divide the domain size
    int nX = static_cast<int>(floor(domainSize[0] / cutoff));
    int nY = static_cast<int>(floor(domainSize[1] / cutoff));
    int nZ = static_cast<int>(floor(domainSize[2] / cutoff));

    // minimum number of cells is 1
    nX = std::max(nX, 1);
    nY = std::max(nY, 1);
    nZ = std::max(nZ, 1);

    // number of cells in each direction
    this->nCells = {nX, nY, nZ};

    // cell size
    double size_x = this->domainSize[0] / nX;
    double size_y = this->domainSize[1] / nY;
    double size_z = this->domainSize[2] / nZ;

    // create cells including halo cells
    for (int i = -1; i <= nZ; i++) {
        for (int j = -1; j <= nY; j++) {
            for (int k = -1; k <= nX; k++) {
                // compute the position of the cell
                std::array<int, 3> index3D = {k, j, i};
                std::array<double, 3> position = {k * size_x, j * size_y, i * size_z};
                std::array<double, 3> size = {size_x, size_y, size_z};
                int index = get1DIndex(index3D);
                std::set<int> neighbours = getNeighborCells(index);
                // push the cell into the cells vector
                cells.emplace_back(position, size, neighbours);
                for (int l = 0; l < 3; ++l) {
                    if (index3D[l] == -1) {
                        haloCells[2 * l].push_back(index);
                    } else if (index3D[l] == nCells[l]) {
                        haloCells[2 * l + 1].push_back(index);
                    }
                }
                if (k >= 0 && k < nCells[0] && j >= 0 && j < nCells[1] && i >= 0 && i < nCells[2]) {
                    domainCells.push_back(index);
                }
            }
        }
    }
    // assign particles to cells

    for(int i = 0; i < particles.size(); i++){

        const auto& pos = particles[i].getX(); // Get particle position
        
        // Add particle index to the corresponding cell
        int idx = getCellIndex(pos);
        std::array<int, 3> index3D = get3DIndex(idx);
        if(index3D[0] >= 0 && index3D[0] < nCells[0] && index3D[1] >= 0 && index3D[1] < nCells[1] && index3D[2] >= 0 && index3D[2] < nCells[2]) {
            cells[idx].addIndex(i);
        } else if (particles[i].isInDomain()) {
            (this->particles)[i].removeFromDomain();
            particleNumber--;
        }
    }
    spdlog::trace("LinkedCellContainer generated!");


}


void LinkedCellContainer::updateV(double delta_t) {

    #pragma omp parallel default(none) shared(particles,cells, delta_t) num_threads(THREADS)
    {
        #pragma omp for schedule(guided, 4)
        for (int i: domainCells) {
            for (int p: cells[i].getParticleIndices()) {
                particles[p].updateV(delta_t);
            }
        }
    }
}

void LinkedCellContainer::updateF() {

    resetF();

    spdlog::trace("Try to use {} threads", THREADS);

    #pragma omp parallel default(none) shared( particles, cells, domainCells, cutoff, force) num_threads(THREADS)
    {
        #pragma omp for schedule(guided, 4)
        for(int i: domainCells){

            std::vector<int> pointCellParticles = cells[i].getParticleIndices();

            // update forces between neighbouring cells
            auto neighbors = cells[i].getNeighbours();

            std::vector<int> neighborList(neighbors.begin(), neighbors.end());
            
            // update forces between neighboring cells
            for(int j = 0 ; j < neighborList.size(); j++){
                if(i < neighborList[j]){
                    updateFCells(i, neighborList[j]);
                }
            }
            

            for(unsigned long j = 0; j < pointCellParticles.size(); ++j){
                for(unsigned long k = j + 1; k < pointCellParticles.size(); ++k){

                    if (particles[pointCellParticles[j]].isStationary() && particles[pointCellParticles[k]].isStationary()) continue;

                    double distSquare = ArrayUtils::L2NormSquare(particles[pointCellParticles[j]].getX() - particles[pointCellParticles[k]].getX());

                    // if the distance is greater than the cutoff, skip the calculation
                    if(distSquare > cutoff * cutoff) continue;

                    std::array<double, 3> forceIJ = force->force(particles[pointCellParticles[j]], particles[pointCellParticles[k]]);
                    

                    particles[pointCellParticles[j]].lock();
                    particles[pointCellParticles[j]].addForce(-1 * forceIJ);
                    particles[pointCellParticles[j]].unlock();

                    particles[pointCellParticles[k]].lock();
                    particles[pointCellParticles[k]].addForce(forceIJ);
                    particles[pointCellParticles[k]].unlock();

                    
                }
            }
            
            

        }
    }

    // for(int thread_id = 0; thread_id < THREADS; thread_id++){
    //     for (size_t i = 0; i < particles.size(); ++i) {
    //         particles[i].addForce(thread_forces[thread_id][i]);
    //     }
    //     thread_forces[thread_id].assign(particles.size(), {0.0, 0.0, 0.0});
    // }
        
}


void LinkedCellContainer::updateFCells(int c1, int c2) {

    std::vector<int> v1 = cells[c1].getParticleIndices();
    std::vector<int> v2 = cells[c2].getParticleIndices();   

    std::array<int, 3> idx3D1 = get3DIndex(c1);
    std::array<int, 3> idx3D2 = get3DIndex(c2);

    // Consider periodic boundary condition in the force calculation.
    std::vector<std::array<double, 3>> offsets;

    for (int z = -1 * (boundaryConditions[4] == PERIODIC); z <= (boundaryConditions[5] == PERIODIC); ++z) {
        for (int y = -1 * (boundaryConditions[2] == PERIODIC); y <= (boundaryConditions[3] == PERIODIC); ++y) {
            for (int x = -1 * (boundaryConditions[0] == PERIODIC); x <= (boundaryConditions[1] == PERIODIC); ++x) {
                if (std::abs(idx3D1[0] + x * nCells[0] - idx3D2[0]) <= 1 &&
                    std::abs(idx3D1[1] + y * nCells[1] - idx3D2[1]) <= 1 &&
                    std::abs(idx3D1[2] + z * nCells[2] - idx3D2[2]) <= 1 ) {
                    offsets.push_back({x * domainSize[0], y * domainSize[1], z * domainSize[2]});
                }
            }
        }
    }

    //#pragma omp parallel shared(v1, v2, particles, cutoff, force, offsets, thread_forces) num_threads(InnerThreads)
    {   
        for (int i : v1) {
            for (int j : v2) {
                if (particles[i].isStationary() && particles[j].isStationary()) continue;

                for (std::array<double, 3> offset: offsets) {

                    double distSquare = ArrayUtils::L2NormSquare(particles[i].getX() + offset - particles[j].getX());
                    // if the distance is greater than the cutoff, skip the calculation
                    if(distSquare <= cutoff * cutoff) {

                        std::array<double, 3> pos = particles[i].getX();

                        particles[i].lock();
                        particles[i].setX(pos + offset); 
                        std::array<double, 3> forceIJ = force->force(particles[i], particles[j]);
                        particles[i].setX(pos);

                        particles[i].addForce(-1 * forceIJ);
                        particles[i].unlock();

                        particles[j].lock();
                        particles[j].addForce(forceIJ);
                        particles[j].unlock();
                    }
                }
            }
        }
    }
    
}

void LinkedCellContainer::updateX(double delta_t){

    for (int i = 0; i < particles.size(); i++) {

        if (particles[i].isInDomain()) {
            int cellidx_before = getCellIndex(particles[i].getX());

            particles[i].updateX(delta_t);

            int cellidx_after = getCellIndex(particles[i].getX());

            if(cellidx_before != cellidx_after){
                if (cellidx_before >= 0 && cellidx_before < cells.size()) {
                    cells[cellidx_before].removeIndex(i);
                }
                if(cellidx_after >= 0 && cellidx_after < cells.size()) {
                    cells[cellidx_after].addIndex(i);
                } else {
                    particles[i].removeFromDomain();
                    particleNumber--;
                }
            }
        }

    }
    updateHalo(delta_t);
}

// only used once in the constructor.
std::set<int> LinkedCellContainer::getNeighborCells(int cellIndex) {
    std::array<int, 3> index3D = get3DIndex(cellIndex);
    std::set<int> neighbors;

    // Halo cells don't need neighbours as they are not considered in the force calculation.
    if (isHaloCell(cellIndex)) {
        return neighbors;
    }
    //
    std::array<int, 6> neighbourBounds{};

    for (int i = 0; i < 3; ++i) {
        neighbourBounds[2 * i] = boundaryConditions[2 * i] == PERIODIC ? index3D[i] - 1 : std::max(0, index3D[i] - 1);
        neighbourBounds[2 * i + 1] = boundaryConditions[2 * i + 1] == PERIODIC ? index3D[i] + 1 : std::min(nCells[i] - 1, index3D[i] + 1);
    }

    // get all neighbors, for periodic boundary condition the cells on the other side are considered as well.
    for (int z = neighbourBounds[4]; z <= neighbourBounds[5]; ++z) {
        for (int y = neighbourBounds[2]; y <= neighbourBounds[3]; ++y) {
            for (int x = neighbourBounds[0]; x <= neighbourBounds[1]; ++x) {
                int neighborIndex = get1DIndex({(x + nCells[0]) % nCells[0], (y + nCells[1]) % nCells[1], (z + nCells[2]) % nCells[2]});
                if (neighborIndex != cellIndex) {
                    neighbors.insert(neighborIndex);
                }
            }
        }
    }
    return neighbors;
}


int LinkedCellContainer::get1DIndex(std::array<int, 3> index3D) {
    return (index3D[2] + 1) * (nCells[0] + 2) * (nCells[1] + 2) + (index3D[1] + 1) * (nCells[0] + 2) + index3D[0] + 1;
}

std::array<int, 3> LinkedCellContainer:: get3DIndex(int cellIndex) {
    int idxZ = cellIndex / ((nCells[0] + 2) * (nCells[1] + 2)) - 1;
    int idxY = (cellIndex - (idxZ + 1) * (nCells[0] + 2) * (nCells[1] + 2)) / (nCells[0] + 2) - 1;
    int idxX = cellIndex - (idxZ + 1) * (nCells[0] + 2) * (nCells[1] + 2) - (idxY + 1) * (nCells[0] + 2) - 1;
    return {idxX, idxY, idxZ};
}


bool LinkedCellContainer::isBoundaryCell(int index) {
    std::array<int, 3> Index3D = get3DIndex(index);
    for (int i = 0; i < 3; ++i) {
        if (Index3D[i] == 0 || Index3D[0] == nCells[0] - 1) {
            return true;
        }
    }
    return false;
}

bool LinkedCellContainer::operator==(const LinkedCellContainer &other) const {
    if (getParticleNumber() != other.getParticleNumber()) {
        return false;
    }
    for (int i = 0, i1 = 0, i2 = 0; i < getParticleNumber(); ++i) {
        while (!getParticles()[i1].isInDomain()) i1++;
        while (!other.getParticles()[i2].isInDomain()) i2++;
        if (!(getParticles()[i1++] == other.getParticles()[i2++])) {
            return false;
        }
    }
    if (!(domainSize == other.domainSize) || cutoff != other.cutoff || !(boundaryConditions == other.boundaryConditions)) {
        return false;
    }
    return true;
}

void LinkedCellContainer::addParticle(const Particle& particle) {
    int cellIndex = getCellIndex(particle.getX());
    if (isDomainCell(cellIndex) && particle.isInDomain()) {
        cells[cellIndex].addIndex(static_cast<int>(particles.size()));
        particles.push_back(particle);
        particleNumber++;
    }
}

void LinkedCellContainer::addCluster(const Cluster &cluster) {
    int size_old = static_cast<int>(particles.size());
    cluster.createParticles(particles);
    for (int i = size_old; i < particles.size(); i++) {
        int cellIndex = getCellIndex(particles[i].getX());
        if (isDomainCell(cellIndex)) {
            cells[cellIndex].addIndex(i);
            particleNumber++;
        } else {
            particles[i].removeFromDomain();
        }
    }
}

bool LinkedCellContainer::isHaloCell(int index) {
    if (index < -1 || index > cells.size()) {
        return false;
    }
    std::array<int, 3> index3D = get3DIndex(index);
    return index3D[0] == -1 || index3D[0] == nCells[0] || index3D[1] == -1 || index3D[1] == nCells[1] || index3D[2] == -1 || index3D[2] == nCells[2];
}

bool LinkedCellContainer::isDomainCell(int index) {
    if (index < 0 || index >= cells.size()) {
        return false;
    }
    std::array<int, 3> index3D = get3DIndex(index);
    return index3D[0] >= 0 && index3D[0] < nCells[0] && index3D[1] >= 0 && index3D[1] < nCells[1] && index3D[2] >= 0 && index3D[2] < nCells[2];
}



void LinkedCellContainer::removeFromHalo(Direction direction) {
    for(int i : haloCells[direction]) {
        for (int p: cells[i].getParticleIndices()) {
            if (particles[p].isInDomain()) {
                particles[p].removeFromDomain();
                particleNumber--;
            }
        }
        cells[i].clear();
    }
}

void LinkedCellContainer::updateHalo(Direction direction, BoundaryCondition boundaryCondition, double deltaT) {
    for(int i : haloCells[direction]) {
        for (int p : cells[i].getParticleIndices()) {
            std::array<double, 3> pos = particles[p].getX();
            std::array<double, 3> vel = particles[p].getV();
            switch (boundaryCondition) {
                case OUTFLOW:
                    if (particles[p].isInDomain()) {
                        particles[p].removeFromDomain();
                        particleNumber--;
                    }
                    break;
                case REFLECTING:

                    // change sign of the corresponding velocity component.
                    // Also flipping the term of f_{old} in the next update of the velocity.
                    // Note: Since the next update will add f_old / (2m) * dt to the velocity, subtracting twice of the
                    // term in advance effectively flips this term

                    vel[direction / 2] *= -1;
                    vel[direction / 2] -= particles[p].getF()[direction / 2] * deltaT / particles[p].getM();
                    particles[p].setV(vel);

                    // reflect the position of the particle with respect to the boundary.

                    pos[direction / 2] = 2 * (direction % 2) * domainSize[direction / 2] - pos[direction / 2];
                    particles[p].setX(pos);
                    cells[getCellIndex(pos)].addIndex(p);
                    break;
                case PERIODIC:

                    // Move the particle to the other side once it left the domain.

                    pos[direction / 2] = pos[direction / 2] + pow(-1, direction % 2) * domainSize[direction / 2];
                    particles[p].setX(pos);
                    cells[getCellIndex(pos)].addIndex(p);
            }
        }
        cells[i].clear();
    }
}

void LinkedCellContainer::updateHalo(double deltaT) {
    for (int i = 0; i < 6; ++i) {
        updateHalo(static_cast<Direction>(i), boundaryConditions[i], deltaT);
    }
}

const std::vector<Cell> &LinkedCellContainer::getCells() const {
    return cells;
}

double LinkedCellContainer::getCutoff() const {
    return cutoff;
}

const std::array<double, 3> &LinkedCellContainer::getDomainSize() const {
    return domainSize;
}

const std::array<int, 3> &LinkedCellContainer::getNCells() const {
    return nCells;
}

const std::array<BoundaryCondition, 6> &LinkedCellContainer::getBoundaryConditions() const {
    return boundaryConditions;
}

std::string LinkedCellContainer::toString() {
    std::stringstream buf;
    buf << ParticleContainer::toString();
    buf << "Domain size: " << domainSize << std::endl;
    buf << "Cutoff radius: " << cutoff << std::endl;
    buf << "Boundary conditions: " << boundaryConditions << std::endl;
    return buf.str();
}


