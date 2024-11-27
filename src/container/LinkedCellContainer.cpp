#include "LinkedCellContainer.h"
#include "spdlog/spdlog.h"
#include <vector>
#include <array>
#include <iostream>
#include <math.h>
#include "utils/ArrayUtils.h"


LinkedCellContainer::LinkedCellContainer(std::vector<Particle> &particles, std::unique_ptr<Force> &f_ptr, std::array<double, 3> domainSize, double cutoff) :
        Container(particles, f_ptr) {

    this->cutoff = cutoff;
    this->domainSize = domainSize;

    //TODO: Initialize the cells
    // compute cell indicies size
    int nX = static_cast<int>(ceil(std::abs(domainSize[0]) / cutoff));
    int nY = static_cast<int>(ceil(std::abs(domainSize[1]) / cutoff));
    int nZ = static_cast<int>(ceil(std::abs(domainSize[2]) / cutoff));

    // minimum number of cells is 1
    nX = std::max(nX, 1);
    nY = std::max(nY, 1);
    nZ = std::max(nZ, 1);

    this->nCells = {nX, nY, nZ};


    // calculate the start position of the domain
    double startX = domainSize[0] < 0 ? domainSize[0] : 0;
    double startY = domainSize[1] < 0 ? domainSize[1] : 0;
    double startZ = domainSize[2] < 0 ? domainSize[2] : 0;

    // create cells
    for (int i = 0; i < nZ; i++) {
        for (int j = 0; j < nY; j++) {
            for (int k = 0; k < nX; k++) {
                // compute the position of the cell
                std::array<double, 3> position = { startX + k * cutoff,  startY + j * cutoff, startZ + i * cutoff };
                std::array<double, 3> size = {cutoff, cutoff, cutoff};
                // push the cell into the cells vector
                cells.emplace_back(Cell(position, size));
            }
        }
    }

    // assign particles to cells    
    for(int i = 0; i< particles.size(); i++){

        const auto& pos = particles[i].getX(); // Get particle position

        // Add particle index to the corresponding cell
        cells[getCellIndex(pos)].addIndex(i);
    } 
    spdlog::trace("LinkedCellContainer generated!");
}


//TODO: Implement the function to update the velocities
void LinkedCellContainer::updateV(double delta_t) {
        for (auto &p : particles) {
        std::array<double, 3> vec = p.getV() + (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
        p.setV(vec);
    }
}

//TODO: Implement the function to update the forces
void LinkedCellContainer::updateF(bool newton3) {

}

void LinkedCellContainer::updateX(double delta_t){
    for (auto & particle : particles) {
        std::array<double, 3> vec = particle.getX() + delta_t * (particle.getV() + (delta_t / (2 * particle.getM())) * particle.getF());
        particle.setX(vec);
  }
}


std::vector<int> LinkedCellContainer::getNeighborCells(int cellIndex) {
    int idxZ = cellIndex / (nCells[0] * nCells[1]);
    int idxY = (cellIndex - idxZ * nCells[0] * nCells[1]) / nCells[0];
    int idxX = cellIndex - idxZ * nCells[0] * nCells[1] - idxY * nCells[0];
    std::vector<int> neighbors;

    for (int i = std::max(0, idxZ - 1); i <= std::min(idxZ + 1, nCells[2] - 1); i++) {
        for (int j = std::max(0, idxY - 1); j <= std::min(idxY + 1, nCells[1] - 1); j++) {
            for (int k = std::max(0, idxX - 1); k <= std::min(idxX + 1, nCells[0] - 1); k++) {

                neighbors.push_back(i * nCells[0] * nCells[1] + j * nCells[0] + k);
            }
        }
    }

    return neighbors;
}


int LinkedCellContainer::getCellIndex(std::array<double, 3> positions) {

    // calculate the start position of the domain
    double startX = domainSize[0] < 0 ? domainSize[0] : 0;
    double startY = domainSize[1] < 0 ? domainSize[1] : 0;
    double startZ = domainSize[2] < 0 ? domainSize[2] : 0;

    // Calculate cell positions for each dimension
    int idxX = static_cast<int>(std::floor((positions[0] - startX) / cutoff));
    int idxY = static_cast<int>(std::floor((positions[1]- startY) / cutoff));
    int idxZ = static_cast<int>(std::floor((positions[2] - startZ) / cutoff));

    // Check if the position is exactly divisible by cutoff and adjust positions
    if (std::fmod(positions[0] - startX, cutoff) == 0 && idxX > 0) idxX -= 1;
    if (std::fmod(positions[1] - startY, cutoff) == 0 && idxY > 0) idxY -= 1;
    if (std::fmod(positions[2] - startZ, cutoff) == 0 && idxZ > 0) idxZ -= 1;

    // Calculate the linear cell index
    int cellIndex = idxZ * nCells[0] * nCells[1] + idxY * nCells[0] + idxX;

    return cellIndex;

}