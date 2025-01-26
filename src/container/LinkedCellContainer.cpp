#include <omp.h>
#include <vector>
#include <array>
#include <spdlog/spdlog.h>
#include "container/LinkedCellContainer.h"
#include "utils/ArrayUtils.h"
#include "cmath"


LinkedCellContainer::LinkedCellContainer(std::vector<Particle>& particles, std::unique_ptr<Force> &f,
                                         std::array<double, 3> domainSize, double cutoff, std::array<BoundaryCondition, 6> boundaryConditions) :
        ParticleContainer(particles, f), cutoff(cutoff), domainSize(domainSize), boundaryConditions(boundaryConditions){

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

    initializePairs();

    spdlog::trace("LinkedCellContainer generated!");
}


void LinkedCellContainer::updateV(double delta_t) {
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) default(none) shared(particles, cells, delta_t)
#endif
    for (int i = 0; i < domainCells.size(); i++) {
        for (int p: cells[domainCells[i]].getParticleIndices()) {
            particles[p].updateV(delta_t);
        }
    }
}

void LinkedCellContainer::updateF(int strategy) {

    resetF();

    // update forces between neighbouring cells
    if (!strategy) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(guided, 4) default(none) shared(particles, cells, domainCells, cutoff, force)
#endif
        for (int i = 0; i < domainCells.size(); i++) {
            auto neighbors = cells[domainCells[i]].getNeighbours();
            for (auto &neighbor: neighbors) {
                // make sure every pair of cells is only considered once
                if (domainCells[i] < neighbor) {
                    updateFCells(neighbor, domainCells[i], true);
                }
            }
        }

    } else {
        for (const std::vector<Pair>& pairs : cellPairs) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic)
#endif
            for (int i = 0; i < pairs.size(); i++) {
                updateFCells(pairs[i].first, pairs[i].second);
            }
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(guided, 4) default(none) shared(particles, cells, domainCells, cutoff, force)
#endif
        // update forces within a cell.
    for (int i = 0; i < domainCells.size(); i++) {
        std::vector<int> pointCellParticles = cells[domainCells[i]].getParticleIndices();

        for (unsigned long j = 0; j < pointCellParticles.size(); ++j) {
            for (unsigned long k = j + 1; k < pointCellParticles.size(); ++k) {

                if (particles[pointCellParticles[j]].isStationary() &&
                    particles[pointCellParticles[k]].isStationary())
                    continue;

                double distSquare = ArrayUtils::L2NormSquare(
                        particles[pointCellParticles[j]].getX() - particles[pointCellParticles[k]].getX());

                // if the distance is greater than the cutoff, skip the calculation
                if (distSquare > cutoff * cutoff) continue;

                std::array<double, 3> forceIJ = force->force(particles[pointCellParticles[j]],
                                                             particles[pointCellParticles[k]]);

                particles[pointCellParticles[j]].addForce(-1 * forceIJ);
                particles[pointCellParticles[k]].addForce(forceIJ);
            }
        }

    }
}

void LinkedCellContainer::updateFCells(int c1, int c2, bool synchronized){

    std::vector<int> v1 = cells[c1].getParticleIndices();
    std::vector<int> v2 = cells[c2].getParticleIndices();

    std::array<int, 3> idx3D1 = get3DIndex(c1);
    std::array<int, 3> idx3D2 = get3DIndex(c2);

    // Consider periodic boundary condition in the force calculation.
    std::vector<std::array<double, 3>> offsets{};
    std::array<std::vector<double>, 3> offsetsDirection{};

    for (int i = 0; i < 3; ++i) {
        if (std::abs(idx3D1[i] - idx3D2[i]) <= 1) {
            offsetsDirection[i].push_back(0);
        }
        if (std::abs(idx3D1[i] - idx3D2[i]) == nCells[i] - 1 && boundaryConditions[2 * i] == PERIODIC) {
            if (idx3D1[i] >= idx3D2[i]) {
                offsetsDirection[i].push_back(-domainSize[i]);
            }
            if (idx3D1[i] <= idx3D2[i]) {
                offsetsDirection[i].push_back(domainSize[i]);
            }
        }
    }

    for (double offsetX: offsetsDirection[0]) {
        for (double offsetY: offsetsDirection[1]) {
            for (double offsetZ: offsetsDirection[2]) {
                offsets.push_back({offsetX, offsetY, offsetZ});
            }
        }
    }

    for (int i : v1) {
        for (int j : v2) {

            if (particles[i].isStationary() && particles[j].isStationary()) continue;

            for (std::array<double, 3> offset: offsets) {

                double distSquare = ArrayUtils::L2NormSquare(particles[i].getX() + offset - particles[j].getX());
                // if the distance is greater than the cutoff, skip the calculation
                if(distSquare <= cutoff * cutoff) {

                    std::array<double, 3> forceIJ{};

                    std::array<double, 3> pos = particles[i].getX();

                    // in most cases, offset is 0. Hence, this logic avoids copying/moving particle if this is the case.
                    if(offset[0] == 0 && offset[1] == 0 && offset[2] == 0) {
                        forceIJ = force->force(particles[i], particles[j]);
                    } else {
                        // If synchronization is required, then copying can't be avoided, since moving particle might lead to race
                        // conditions, or otherwise the computation of the force must be executed exclusively as well.
                        // In case no synchronization is needed, moving instead of copying causes less overhead.
                        if (synchronized) {
                            Particle mock(particles[i]);
                            mock.setX(pos + offset);
                            forceIJ = force->force(mock, particles[j]);
                        } else {
                            particles[i].setX(pos + offset);
                            forceIJ = force->force(particles[i], particles[j]);
                            particles[i].setX(pos);
                        }
                    }
#ifdef _OPENMP
                    if (synchronized) {
                        particles[i].lock();
                    }
                    particles[i].addForce( -1 * forceIJ);
                    if (synchronized) {
                        particles[i].unlock();
                        particles[j].lock();
                    }
                    particles[j].addForce(forceIJ);
                    if (synchronized) {
                        particles[j].unlock();
                    }
#else
                    particles[i].addForce( -1 * forceIJ);
                    particles[j].addForce(forceIJ);
#endif
                }
            }
        }
    }
}

void LinkedCellContainer::updateX(double delta_t){
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
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
#ifdef _OPENMP
                    #pragma omp critical
#endif
                    {
                        particleNumber--;
                    }

                }
            }
        }

    }
    updateHalo(delta_t);
}

void LinkedCellContainer::initializePairs() {
    for (int deltaZ = -1; deltaZ <= 1; ++deltaZ) {
        for (int deltaY = -1; deltaY <= 1; ++deltaY) {
            std::vector<Pair> pairsRight1{};
            std::vector<Pair> pairsRight2{};

            for (int x = 0; x < nCells[0] - 1; x += 2) {
                for (int y = 0; y < nCells[1]; ++y) {
                    for (int z = 0; z < nCells[2]; ++z) {
                        int first1 = get1DIndex({x, y, z});
                        int second1 = get1DIndex(moduloCellNumber({x + 1, y + deltaY, z + deltaZ}));
                        int first2 = get1DIndex({x + 1, y, z});
                        int second2 = get1DIndex(moduloCellNumber({x + 2, y + deltaY, z + deltaZ}));
                        if (isDomainCell(second1)) {
                            pairsRight1.push_back(Pair{first1, second1});
                        }
                        if (isDomainCell(second2)) {
                            pairsRight2.push_back(Pair{first2, second2});
                        }
                    }
                }
            }
            cellPairs.push_back(pairsRight1);
            cellPairs.push_back(pairsRight2);

            if (boundaryConditions[0] == PERIODIC && nCells[0] % 2 == 1) {
                std::vector<Pair> pairsBoundary{};
                for (int y = 0; y < nCells[1]; ++y) {
                    for (int z = 0; z < nCells[2]; ++z) {
                        int first = get1DIndex({nCells[0] - 1, y, z});
                        int second = get1DIndex(moduloCellNumber({0, y + deltaY, z + deltaZ}));
                        if (isDomainCell(second)) {
                            pairsBoundary.push_back(Pair{first, second});
                        }
                    }
                }
                cellPairs.push_back(pairsBoundary);
            }
        }

        std::vector<Pair> pairsUp1;
        std::vector<Pair> pairsUp2;

        for (int y = 0; y < nCells[1] - 1; y += 2) {
            for (int z = 0; z < nCells[2]; ++z) {
                for (int x = 0; x < nCells[0]; ++x) {
                    int first1 = get1DIndex({x, y, z});
                    int second1 = get1DIndex(moduloCellNumber({x, y + 1, z + deltaZ}));
                    int first2 = get1DIndex({x, y + 1, z});
                    int second2 = get1DIndex(moduloCellNumber({x, y + 2, z + deltaZ}));
                    if (isDomainCell(second1)) {
                        pairsUp1.push_back(Pair{first1, second1});
                    }
                    if (isDomainCell(second2)) {
                        pairsUp2.push_back(Pair{first2, second2});
                    }
                }
            }
        }

        cellPairs.push_back(pairsUp1);
        cellPairs.push_back(pairsUp2);

        if (boundaryConditions[2] == PERIODIC && nCells[1] % 2 == 1) {
            std::vector<Pair> pairsBoundary{};
            for (int x = 0; x < nCells[0]; ++x) {
                for (int z = 0; z < nCells[2]; ++z) {
                    int first = get1DIndex({x, nCells[1] - 1, z});
                    int second = get1DIndex(moduloCellNumber({x, 0, z + deltaZ}));
                    if (isDomainCell(second)) {
                        pairsBoundary.push_back(Pair{first, second});
                    }
                }
            }
            cellPairs.push_back(pairsBoundary);
        }
    }
    std::vector<Pair> pairsFront1{};
    std::vector<Pair> pairsFront2{};

    for (int z = 0; z < nCells[2] - 1; z += 2) {
        for (int x = 0; x < nCells[0]; ++x) {
            for (int y = 0; y < nCells[1]; ++y) {
                int first1 = get1DIndex({x, y, z});
                int second1 = get1DIndex(moduloCellNumber({x, y, z + 1}));
                int first2 = get1DIndex({x, y, z + 1});
                int second2 = get1DIndex(moduloCellNumber({x, y, z + 2}));
                if (isDomainCell(second1)) {
                    pairsFront1.push_back(Pair{first1, second1});
                }
                if (isDomainCell(second2)) {
                    pairsFront2.push_back(Pair{first2, second2});
                }
            }
        }
    }

    cellPairs.push_back(pairsFront1);
    cellPairs.push_back(pairsFront2);

    if (boundaryConditions[4] == PERIODIC && nCells[2] % 2 == 1) {
        std::vector<Pair> pairsBoundary{};
        for (int x = 0; x < nCells[0]; ++x) {
            for (int y = 0; y < nCells[1]; ++y) {
                int first = get1DIndex({x, y, nCells[2] - 1});
                int second = get1DIndex(moduloCellNumber({x, y, 0}));
                if (isDomainCell(second)) {
                    pairsBoundary.push_back(Pair{first, second});
                }
            }
        }
        cellPairs.push_back(pairsBoundary);
    }
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
    if (index3D[0] < -1 || index3D[0] > nCells[0] || index3D[1] < -1 || index3D[1] > nCells[1] || index3D[2] < -1 || index3D[2] > nCells[2]) {
        return -1;
    }
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
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for(int i = 0; i < haloCells[direction].size(); i++) {
        for (int p : cells[haloCells[direction][i]].getParticleIndices()) {
            std::array<double, 3> pos = particles[p].getX();
            std::array<double, 3> vel = particles[p].getV();
            switch (boundaryCondition) {
                case OUTFLOW:
                    if (particles[p].isInDomain()) {
                        particles[p].removeFromDomain();
#ifdef _OPENMP
                        #pragma omp critical
#endif
                        {
                            particleNumber--;
                        }
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
        cells[haloCells[direction][i]].clear();
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

const std::vector<std::vector<Pair>> &LinkedCellContainer::getCellPairs() const {
    return cellPairs;
}

std::string LinkedCellContainer::toString() {
    std::stringstream buf;
    buf << ParticleContainer::toString();
    buf << "Domain size: " << domainSize << std::endl;
    buf << "Cutoff radius: " << cutoff << std::endl;
    buf << "Boundary conditions: " << boundaryConditions << std::endl;
    return buf.str();
}


