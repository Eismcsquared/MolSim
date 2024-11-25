#include "LinkedCellContainer.h"

//TODO: Implement the function to update the velocities
void LinkedCellContainer::updateV(double delta_t) {

}

//TODO: Implement the function to update the forces
void LinkedCellContainer::updateF() {

}

//TODO: Implement the function to update the positions
void LinkedCellContainer::updateX(double delta_t) {

}

//TODO: Implement the function to return the number of particles
unsigned long LinkedCellContainer::getParticleNumber() const {
    return 0;
}

std::unique_ptr<Iterator> LinkedCellContainer::iterator() const {
    return std::make_unique<LinkedCellContainerIterator>(cells.begin(), cells.end());
}

LinkedCellContainerIterator::LinkedCellContainerIterator(std::vector<Cell>::iterator currentCell,
                                                         std::vector<Cell>::iterator end):
                                                         currentCell(currentCell), end(end){
    while ((*currentCell).getParticles().size() == 0) {
        ++currentCell;
    }
    if (currentCell != end) {
        std::vector<Particle> p = (*currentCell).getParticles();
        currentParticle = p.begin();
        endCurrentCell = p.end();
    }
}

Particle &LinkedCellContainerIterator::next() {
    Particle& current = *currentParticle;
    ++currentParticle;
    if (currentParticle == endCurrentCell) {
        ++currentCell;
        while ((*currentCell).getParticles().size() == 0) {
            ++currentCell;
        }
        if (hasNext()) {
            std::vector<Particle> p = (*currentCell).getParticles();
            currentParticle = p.begin();
            endCurrentCell = p.end();
        }
    }
    return current;
}

bool LinkedCellContainerIterator::hasNext() {
    return currentCell != end;
}
