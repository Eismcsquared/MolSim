#include "container/LinkedCellContainer.h"
#include "spdlog/spdlog.h"
#include <vector>
#include <array>
#include <limits>
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

    // create cells
    for (int i = 0; i < nZ; i++) {
        for (int j = 0; j < nY; j++) {
            for (int k = 0; k < nX; k++) {
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
        unsigned int idx = getCellIndex(pos);
        if(idx < cells.size()) {
            cells[idx].addIndex(i);
        }
    }
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

    for(int i =0 ; i< cells.size(); ++i){
        std::vector<unsigned int> pointCellParticles = cells[i].getParticleIndices();

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

void LinkedCellContainer::updateCellF(const std::vector<unsigned int> &v1, const std::vector<unsigned int> &v2, bool newton3){
    for (unsigned long i = 0; i < v1.size(); ++i) {
        for (unsigned long j = 0; j < v2.size(); ++j) {

            //////////////////////////////////////////////////////////
            double dist = ArrayUtils::L2Norm((*particles)[v1[i]].getX() - (*particles)[v2[j]].getX());
            if(dist > cutoff) continue; 
            // if the distance is greater than the cutoff, skip the calculation
            
            if (newton3) {
                std::array<double, 3> forceIJ = f->force((*particles)[v1[i]], (*particles)[v2[j]]);

                std::array<double, 3> v1_force;
                std::array<double, 3> v2_force;

                for(int k = 0; k < 3; k++){
                    v1_force[k] = (*particles)[v1[i]].getF()[k] - forceIJ[k];
                    v2_force[k] = (*particles)[v2[j]].getF()[k] + forceIJ[k];
                }

                (*particles)[v1[i]].setF(v1_force);
                (*particles)[v2[j]].setF(v2_force);
            } else {
                std::array<double, 3> forceIJ = f->force((*particles)[v1[i]], (*particles)[v2[j]]);
                std::array<double, 3> forceJI = f->force((*particles)[v2[j]], (*particles)[v1[i]]);

                std::array<double, 3> v1_force;
                std::array<double, 3> v2_force;

                for(int k = 0; k < 3; k++){
                    v1_force[k] = (*particles)[v1[i]].getF()[k] + forceJI[k];
                    v2_force[k] = (*particles)[v2[j]].getF()[k] + forceIJ[k];
                }

                (*particles)[v1[i]].setF(v1_force);
                (*particles)[v2[j]].setF(v2_force);
            }

        }
    }
}

void LinkedCellContainer::updateX(double delta_t){

    for (int i = 0; i < particles->size(); i++) {

        unsigned int cellidx_before = getCellIndex((*particles)[i].getX());

        std::array<double, 3> vec = (*particles)[i].getX() + delta_t * ((*particles)[i].getV() + (delta_t / (2 * (*particles)[i].getM())) * (*particles)[i].getF());
        (*particles)[i].setX(vec);

        unsigned int cellidx_after = getCellIndex((*particles)[i].getX());


        if(cellidx_before != cellidx_after){
            cells[cellidx_before].removeIndex(i);
            if(cellidx_after >= 0 && cellidx_after < cells.size()) {
                cells[cellidx_after].addIndex(i);
            } 
            else {
                if(boundaryConditions[0] == BoundaryCondition::PERIODIC){
                    (*particles)[i].removeFromDomain();
                } 
                else {
                    
                }
            }
        }
    }
}


std::vector<int> LinkedCellContainer::getNeighborCells(int cellIndex) {
    int idxZ = cellIndex / (nCells[0] * nCells[1]);
    int idxY = (cellIndex - idxZ * nCells[0] * nCells[1]) / nCells[0];
    int idxX = cellIndex - idxZ * nCells[0] * nCells[1] - idxY * nCells[0];
    std::vector<int> neighbors;

    // get all neighbors
    for (int z = std::max(0, idxZ - 1); z <= std::min(idxZ + 1, nCells[2] - 1); ++z) {
        for (int y = std::max(0, idxY - 1); y <= std::min(idxY + 1, nCells[1] - 1); ++y) {
            for (int x = std::max(0, idxX - 1); x <= std::min(idxX + 1, nCells[0] - 1); ++x) {
                int neighborIndex = z * nCells[0] * nCells[1] + y * nCells[0] + x;
                if (neighborIndex != cellIndex) {
                    neighbors.push_back(neighborIndex);
                }
            }
        }
    }

    return neighbors;
}


unsigned int LinkedCellContainer::getCellIndex(std::array<double, 3> positions) {

    // Check if the position is within the domain
    if(positions[0] < 0 || positions[0] > domainSize[0] || positions[1] < 0 ||
    positions[1] > domainSize[1] || positions[2] < 0 || positions[2] > domainSize[2]){
        return std::numeric_limits<unsigned int>::max();
    }

    int idxX = static_cast<int>(positions[0] / cells[0].getSize()[0]);
    int idxY = static_cast<int>(positions[1] / cells[0].getSize()[1]);
    int idxZ = static_cast<int>(positions[2] / cells[0].getSize()[2]);

    // Check if the position is exactly divisible by cutoff and adjust positions
    if (positions[0] == idxX * cells[0].getSize()[0] && idxX > 0) idxX -= 1;
    if (positions[1] == idxY * cells[0].getSize()[1] && idxY > 0) idxY -= 1;
    if (positions[2] == idxZ * cells[0].getSize()[2] && idxZ > 0) idxZ -= 1;


    // Calculate the linear cell index
    int cellIndex = idxZ * nCells[0] * nCells[1] + idxY * nCells[0] + idxX;

    return cellIndex;

}

std::array<unsigned int, 3> LinkedCellContainer:: get3DIndex (unsigned int cellIndex) {
    unsigned int idxZ = cellIndex / (nCells[0] * nCells[1]);
    unsigned int idxY = (cellIndex - idxZ * nCells[0] * nCells[1]) / nCells[0];
    unsigned int idxX = cellIndex - idxZ * nCells[0] * nCells[1] - idxY * nCells[0];
    return {idxX, idxY, idxZ};
}


bool LinkedCellContainer::isBoundaryCell(unsigned int index) {
    std::array<unsigned int, 3> Index3D = get3DIndex(index);
    for (int i = 0; i < 3; ++i) {
        if (Index3D[i] == 0 || Index3D[0] == nCells[0] - 1) {
            return true;
        }
    }
    return false;
}
