#include "container/LinkedCellContainer.h"
#include "spdlog/spdlog.h"
#include <vector>
#include <array>
#include "utils/ArrayUtils.h"
#include "cmath"


LinkedCellContainer::LinkedCellContainer(std::unique_ptr<std::vector<Particle>>& particles, std::unique_ptr<Force>& f,
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
                std::array<double, 3> position = {k * size_x, j * size_y, i * size_z};
                std::array<double, 3> size = {size_x, size_y, size_z};
                // push the cell into the cells vector
                cells.emplace_back(position, size);
            }
        }
    }
    // assign particles to cells

    for(int i = 0; i < this->particles->size(); i++){

        const auto& pos = (*(this->particles))[i].getX(); // Get particle position
        
        // Add particle index to the corresponding cell
        int idx = getCellIndex(pos);
        if(idx >= 0 && idx < cells.size()) {
            cells[idx].addIndex(i);
        }
    }
    ParticleContainer::updateF();
    spdlog::trace("LinkedCellContainer generated!");
}


void LinkedCellContainer::updateV(double delta_t) {
    for (auto &p : *particles) {
        std::array<double, 3> vec = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
        p.setV(vec);
    }
}

void LinkedCellContainer::updateF(bool newton3) {

    for(auto & p1 : *particles){
        p1.setOldF(p1.getF());
        p1.setF({0, 0, 0}); // reset the force
    }

    for(int i =0 ; i < cells.size(); ++i){
        std::vector<int> pointCellParticles = cells[i].getParticleIndices();

        // update forces between neighbouring cells
        auto neighbors = getNeighborCells(i);
        for(auto & neighbor : neighbors){
            // make sure every pair of cells is only considered once
            if (i < neighbor) {
                updateCellF(pointCellParticles, cells[neighbor].getParticleIndices(), newton3);
            }
        }

        // update forces within a cell.
        for(unsigned long j = 0; j < pointCellParticles.size(); ++j){
            for(unsigned long k = j + 1; k < pointCellParticles.size(); ++k){

                //////////////////////////////////////////////////////////
                double dist = ArrayUtils::L2Norm((*particles)[pointCellParticles[j]].getX() - (*particles)[pointCellParticles[k]].getX());

                // if the distance is greater than the cutoff, skip the calculation
                if(dist > cutoff) continue;

                
                std::array<double, 3> forceIJ = f->force((*particles)[pointCellParticles[j]], (*particles)[pointCellParticles[k]]);
                
                std::array<double, 3> v1_force;
                std::array<double, 3> v2_force;
                
                for(int l = 0; l < 3; l++){
                    v1_force[l] = (*particles)[pointCellParticles[j]].getF()[l] - forceIJ[l];
                    v2_force[l] = (*particles)[pointCellParticles[k]].getF()[l] + forceIJ[l];
                }

                (*particles)[pointCellParticles[j]].setF(v1_force);
                (*particles)[pointCellParticles[k]].setF(v2_force);
            }
        }

        
    }

}

void LinkedCellContainer::updateCellF(const std::vector<int> &v1, const std::vector<int> &v2, bool newton3){
    for (int i : v1) {
        for (int j : v2) {

            //////////////////////////////////////////////////////////
            double dist = ArrayUtils::L2Norm((*particles)[i].getX() - (*particles)[j].getX());
            if(dist > cutoff) continue; 
            // if the distance is greater than the cutoff, skip the calculation
            
            if (newton3) {
                std::array<double, 3> forceIJ = f->force((*particles)[i], (*particles)[j]);

                std::array<double, 3> v1_force;
                std::array<double, 3> v2_force;

                for(int k = 0; k < 3; k++){
                    v1_force[k] = (*particles)[i].getF()[k] - forceIJ[k];
                    v2_force[k] = (*particles)[j].getF()[k] + forceIJ[k];
                }

                (*particles)[i].setF(v1_force);
                (*particles)[j].setF(v2_force);
            } else {
                std::array<double, 3> forceIJ = f->force((*particles)[i], (*particles)[j]);
                std::array<double, 3> forceJI = f->force((*particles)[j], (*particles)[i]);

                std::array<double, 3> v1_force;
                std::array<double, 3> v2_force;

                for(int k = 0; k < 3; k++){
                    v1_force[k] = (*particles)[i].getF()[k] + forceJI[k];
                    v2_force[k] = (*particles)[j].getF()[k] + forceIJ[k];
                }

                (*particles)[i].setF(v1_force);
                (*particles)[j].setF(v2_force);
            }

        }
    }
}

void LinkedCellContainer::updateX(double delta_t){

    for (int i = 0; i < particles->size(); i++) {
        int cellidx_before = getCellIndex((*particles)[i].getX());

        std::array<double, 3> vec = (*particles)[i].getX() + delta_t * ((*particles)[i].getV() + (delta_t / (2 * (*particles)[i].getM())) * (*particles)[i].getF());
        (*particles)[i].setX(vec);

        int cellidx_after = getCellIndex((*particles)[i].getX());


        if(cellidx_before != cellidx_after){
            cells[cellidx_before].removeIndex(i);
            if(cellidx_after >= 0 && cellidx_after < cells.size()) {
                cells[cellidx_after].addIndex(i);
            }
        }
    }
    updateHalo();
}


std::vector<int> LinkedCellContainer::getNeighborCells(int cellIndex) {
    std::array<int, 3> index3D = get3DIndex(cellIndex);
    std::vector<int> neighbors;

    // get all neighbors, also halo cells
    for (int z = std::max(-1, index3D[2] - 1); z <= std::min(index3D[2] + 1, nCells[2]); ++z) {
        for (int y = std::max(-1, index3D[1] - 1); y <= std::min(index3D[1] + 1, nCells[1]); ++y) {
            for (int x = std::max(-1, index3D[0] - 1); x <= std::min(index3D[0] + 1, nCells[0]); ++x) {
                int neighborIndex = get1DIndex({x, y, z});
                if (neighborIndex != cellIndex) {
                    neighbors.push_back(neighborIndex);
                }
            }
        }
    }
    return neighbors;
}


int LinkedCellContainer::getCellIndex(std::array<double, 3> positions) {

    std::array<double, 3> cellSize = cells[0].getSize();
    // Check if the position is within the domain + halo region
    if(positions[0] < -cellSize[0] || positions[0] > domainSize[0] + cellSize[0] || positions[1] < -cellSize[1] ||
    positions[1] > domainSize[1] + cellSize[1] || positions[2] < -cellSize[2] || positions[2] > domainSize[2] + cellSize[2]){
        return -1;
    }

    int idxX = std::floor(positions[0] / cellSize[0]);
    int idxY = std::floor(positions[1] / cellSize[1]);
    int idxZ = std::floor(positions[2] / cellSize[2]);

    // Calculate the linear cell index
    int cellIndex = get1DIndex({idxX, idxY, idxZ});

    return cellIndex;
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
    for (int i = 0; i < getParticleNumber(); ++i) {
        if (!(getParticles()[i] == other.getParticles()[i])) {
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
    if (cellIndex >= 0 && cellIndex < cells.size()) {
        cells[cellIndex].addIndex((int) particles->size());
    }
    particles->push_back(particle);
    updateF(true);
}

void LinkedCellContainer::addCluster(const Cluster &cluster) {
    int size_old = (int) particles->size();
    cluster.createParticles(*particles);
    for (int i = size_old; i < particles->size(); i++) {
        int cellIndex = getCellIndex((*particles)[i].getX());
        if (cellIndex >= 0 && cellIndex < cells.size()) {
            cells[cellIndex].addIndex(i);
        }
    }
    updateF(true);
}

bool LinkedCellContainer::isHaloCell(int index) {
    std::array<int, 3> index3D = get3DIndex(index);
    return index3D[0] < 0 || index3D[0] >= nCells[0] || index3D[1] < 0 || index3D[1] >= nCells[1] || index3D[2] < 0 || index3D[2] >= nCells[2];
}

std::vector<int> LinkedCellContainer::getAllHaloIndices(Direction direction) {
    std::vector<int> haloIndices;
    switch (direction) {
        case LEFT:
            for (int z = -1; z <= nCells[2] ; ++z) {
                for (int y = -1; y <= nCells[1]; ++y) {
                    haloIndices.push_back(get1DIndex({-1, y, z}));
                }
            }
            break;
        case RIGHT:
            for (int z = -1; z <= nCells[2] ; ++z) {
                for (int y = -1; y <= nCells[1]; ++y) {
                    haloIndices.push_back(get1DIndex({nCells[0], y, z}));
                }
            }
            break;
        case DOWN:
            for (int z = -1; z <= nCells[2] ; ++z) {
                for (int x = -1; x <= nCells[0]; ++x) {
                    haloIndices.push_back(get1DIndex({x, -1, z}));
                }
            }
            break;
        case UP:
            for (int z = -1; z <= nCells[2] ; ++z) {
                for (int x = -1; x <= nCells[0]; ++x) {
                    haloIndices.push_back(get1DIndex({x, nCells[1], z}));
                }
            }
            break;
        case BACK:
            for (int y = -1; y <= nCells[1] ; ++y) {
                for (int x = -1; x <= nCells[0]; ++x) {
                    haloIndices.push_back(get1DIndex({x, y, -1}));
                }
            }
            break;
        case FRONT:
            for (int y = -1; y <= nCells[1] ; ++y) {
                for (int x = -1; x <= nCells[0]; ++x) {
                    haloIndices.push_back(get1DIndex({x, y, nCells[2]}));
                }
            }
    }
    return haloIndices;
}

void LinkedCellContainer::removeFromHalo(Direction direction) {
    for(int i : getAllHaloIndices(direction)) {
        for (int p: cells[i].getParticleIndices()) {
            (*particles)[p].removeFromDomain();
            cells[i].removeIndex(p);
        }
    }
}

void LinkedCellContainer::updateHalo(Direction direction, BoundaryCondition boundaryCondition) {
    for(int i : getAllHaloIndices(direction)) {
        for (int p : cells[i].getParticleIndices()) {
            switch (boundaryCondition) {
                case OUTFLOW:
                    (*particles)[p].removeFromDomain();
                    cells[i].removeIndex(p);
                    break;
                case REFLECTING:
                    // change sign of the corresponding velocity component if the particle is flying away.
                    if ((*particles)[p].getV()[direction / 2] * (direction % 2 - 0.5) > 0) {
                        std::array<double, 3> vel = (*particles)[p].getV();
                        vel[direction / 2] *= -1;
                        (*particles)[p].setV(vel);
                    }
                    break;
                case PERIODIC:
                    cells[i].removeIndex(p);
                    std::array<double, 3> pos = (*particles)[p].getX();
                    pos[direction / 2] = pos[direction / 2] + pow(-1, direction % 2) * domainSize[direction / 2];
                    (*particles)[p].setX(pos);
            }
        }
    }
}

void LinkedCellContainer::updateHalo() {
    for (int i = 0; i < 6; ++i) {
        updateHalo(static_cast<Direction>(i), boundaryConditions[i]);
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


