#include "LinkedCellContainer.h"

LinkedCellContainer::LinkedCellContainer(std::vector<Particle> &particles, std::unique_ptr<Force> &f_ptr) :
        Container(particles, f_ptr) {
    //TODO: Initialize the cells
}

//TODO: Implement the function to update the velocities
void LinkedCellContainer::updateV(double delta_t) {

}

//TODO: Implement the function to update the forces
void LinkedCellContainer::updateF(bool newton3) {

}

//TODO: Implement the function to update the positions
void LinkedCellContainer::updateX(double delta_t) {

}

